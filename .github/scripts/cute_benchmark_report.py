#!/usr/bin/env python3
"""Report of a CUTE benchmark run of Uno, optionally compared against another PR.

The results of the current run are read from <results-dir>/results.tsv and
<results-dir>/meta.json (written by the CUTE benchmark workflow). Each row of
results.tsv is
    <name> <solved|failed> <wall seconds> <model read seconds> <solve seconds> <objective evaluations> <iterations>
the last four fields being empty when unknown (no optimization summary in the log).

Without --compare-pr, the report gives the number of solved instances and the
shifted geometric means of the objective evaluations and of the solve time.

With --compare-pr, the "before" results are those of the latest
`<artifact-prefix><compare-pr>` artifact (another run of the same preset). The
report then compares both runs on their common instances and a Dolan-Moré
performance profile with respect to the number of objective evaluations is
written to --profile. The report embeds it as ![...](PROFILE_URL_PLACEHOLDER),
which the workflow substitutes with the URL of the published image. If the
"before" results cannot be found, a warning is printed and the plain report is
written instead.

"comparison=true|false" and "profile=true|false" (whether a profile was drawn) are
appended to $GITHUB_OUTPUT (if set).

Standard library only, except matplotlib for the profile (imported only then).
"""

import argparse
import io
import json
import math
import os
import sys
import zipfile
from typing import NamedTuple, Optional

# reuse the GitHub API helpers of the other benchmarks
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from benchmark_compare import github_request, download_artifact_zip  # noqa: E402

PROFILE_URL_PLACEHOLDER = "PROFILE_URL_PLACEHOLDER"
# GitHub rejects comments longer than 65536 characters
MAX_COMMENT_LENGTH = 60000
# shifts of the shifted geometric means
EVALUATIONS_SHIFT = 10.0
# most CUTE instances solve in milliseconds: a 10s shift (Mittelmann) would average out any change
TIME_SHIFT = 1.0
# names listed in full in the comment (newly solved/failed, only in one run)
MAX_LISTED_NAMES = 200
# relative change in objective evaluations (after/before) beyond which an instance is flagged
EVALUATIONS_RATIO_THRESHOLD = 1.25


class Result(NamedTuple):
    solved: bool
    wall: float
    read: Optional[float]
    solve: Optional[float]
    objective_evaluations: Optional[int]
    iterations: Optional[int]


def optional(field, convert):
    field = field.strip()
    return convert(field) if field else None


def parse_results_tsv(text):
    results = {}
    for line in text.splitlines():
        # not strip(): it would drop the trailing empty fields
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) != 7:
            continue
        name, status, wall, read, solve, evaluations, iterations = fields
        results[name.strip()] = Result(status.strip() == "solved", float(wall), optional(read, float),
                                       optional(solve, float), optional(evaluations, int), optional(iterations, int))
    return results


def fetch_from_artifact(repository, artifact_name, token, current_run_id):
    listing = github_request(f"/repos/{repository}/actions/artifacts?name={artifact_name}&per_page=100", token)
    artifacts = [artifact for artifact in listing.get("artifacts", [])
                 if not artifact.get("expired")
                 and str((artifact.get("workflow_run") or {}).get("id")) != str(current_run_id)]
    if not artifacts:
        return None
    latest = max(artifacts, key=lambda artifact: artifact["created_at"])
    with zipfile.ZipFile(io.BytesIO(download_artifact_zip(repository, latest["id"], token))) as archive:
        results = parse_results_tsv(archive.read("results.tsv").decode())
        meta = json.loads(archive.read("meta.json").decode()) if "meta.json" in archive.namelist() else {}
    return (results, meta) if results else None


def shifted_geometric_mean(values, shift):
    if not values:
        return float("nan")
    return math.exp(sum(math.log(value + shift) for value in values) / len(values)) - shift


def format_number(value, unit="", digits=2):
    return "—" if value is None or math.isnan(value) else f"{value:.{digits}f}{unit}"


def format_time(value):
    return format_number(value, "s", 3)


def format_names(names):
    """Comma-separated names, capped at MAX_LISTED_NAMES to keep the comment postable."""
    listed = ", ".join(f"`{name}`" for name in names[:MAX_LISTED_NAMES])
    return listed + (f", … {len(names) - MAX_LISTED_NAMES} more" if len(names) > MAX_LISTED_NAMES else "")


def annotation_title(title):
    """Workflow-command properties are split on ',' and ':'."""
    return title.replace("%", "%25").replace(",", "%2C").replace(":", "%3A")


def format_ratio(before, after):
    if before is None or after is None or math.isnan(before) or math.isnan(after) or before <= 0:
        return "—"
    return f"{after / before:.2f}×"


def has_counts(result):
    return result.solved and result.objective_evaluations is not None


def run_description(meta):
    parts = []
    if meta.get("pr"):
        parts.append(f"PR #{meta['pr']}")
    if meta.get("sha"):
        parts.append(f"commit `{meta['sha'][:7]}`")
    if meta.get("run_url"):
        parts.append(f"[run]({meta['run_url']})")
    return " · ".join(parts)


def configuration(meta):
    return (f"preset `{meta.get('preset')}` · linear solver `{meta.get('linear_solver')}` · "
            f"time_limit={meta.get('time_limit')}s · {meta.get('selection')}")


def build_summary(title, results, meta, with_table=False):
    """Plain report of a single run (with_table: one row per instance, for the run summary)."""
    solved = [name for name, result in results.items() if result.solved]
    counted = [name for name in solved if has_counts(results[name])]
    timed = [name for name in solved if results[name].solve is not None]
    sgm_evaluations = shifted_geometric_mean([results[name].objective_evaluations for name in counted],
                                             EVALUATIONS_SHIFT)
    sgm_time = shifted_geometric_mean([results[name].solve for name in timed], TIME_SHIFT)
    failed = sorted((name for name, result in results.items() if not result.solved), key=str.casefold)

    lines = [
        f"## {title}",
        "",
        f"{run_description(meta)} · {configuration(meta)}",
        "",
        f"**Solved {len(solved)} / {len(results)} instances**",
        "",
        "| shifted geometric mean (solved instances) | |",
        "|:--|--:|",
        f"| objective evaluations ({len(counted)}) | {format_number(sgm_evaluations)} |",
        f"| solve time ({len(timed)}) | {format_time(sgm_time)} |",
        "",
        f"<sub>shifts: {EVALUATIONS_SHIFT:g} evaluations, {TIME_SHIFT:g}s</sub>",
        "",
    ]
    if with_table:
        lines += [f"<details><summary>Per-instance results ({len(results)} instances)</summary>", "",
                  "| | instance | objective evaluations | iterations | solve time | wall time |",
                  "|:-:|:--|--:|--:|--:|--:|"]
        for name in sorted(results, key=lambda name: (name.casefold(), name)):
            result = results[name]
            lines.append(f"| {'✅' if result.solved else '❌'} | `{name}` | "
                         f"{'—' if result.objective_evaluations is None else result.objective_evaluations} | "
                         f"{'—' if result.iterations is None else result.iterations} | "
                         f"{format_time(result.solve)} | {format_time(result.wall)} |")
        lines += ["", "</details>", ""]
    elif failed:
        lines += [f"<details><summary>Failed instances ({len(failed)})</summary>", "",
                  format_names(failed), "", "</details>", ""]
    console = [f"Solved {len(solved)} / {len(results)} instances",
               f"Shifted geometric mean of the objective evaluations ({len(counted)} solved instances, "
               f"shift {EVALUATIONS_SHIFT:g}): {format_number(sgm_evaluations)}",
               f"Shifted geometric mean of the solve time ({len(timed)} solved instances, "
               f"shift {TIME_SHIFT:g}s): {format_time(sgm_time)}"]
    return "\n".join(lines), console


def performance_profile(before, after, names):
    """Dolan-Moré profile on the objective evaluations: (taus, rho_before, rho_after).
    A failed instance has an infinite ratio; an instance failed by both is never reached."""
    def cost(result):
        return max(result.objective_evaluations, 1) if has_counts(result) else math.inf

    ratios_before, ratios_after = [], []
    for name in names:
        b, a = cost(before[name]), cost(after[name])
        best = min(b, a)
        ratios_before.append(b / best if best < math.inf else math.inf)
        ratios_after.append(a / best if best < math.inf else math.inf)
    finite = [ratio for ratio in ratios_before + ratios_after if ratio < math.inf]
    tau_max = max(2.0, 1.1 * max(finite, default=1.0))
    taus = sorted(set([1.0, tau_max] + finite))
    def rho(ratios):
        return [sum(1 for ratio in ratios if ratio <= tau) / len(ratios) for tau in taus]
    return taus, rho(ratios_before), rho(ratios_after)


def plot_profile(path, taus, rho_before, rho_after, label_before, label_after, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    figure, axes = plt.subplots(figsize=(8, 5), dpi=120)
    axes.step(taus, rho_before, where="post", label=label_before, color="#4269d0", linewidth=2)
    axes.step(taus, rho_after, where="post", label=label_after, color="#ef8a0c", linewidth=2, linestyle="--")
    axes.set_xscale("log", base=2)
    axes.set_xlim(1.0, taus[-1])
    axes.xaxis.set_major_formatter(FuncFormatter(lambda tau, _: f"{tau:g}"))
    axes.set_ylim(0.0, 1.02)
    axes.set_xlabel("τ (ratio to the fewest objective evaluations)")
    axes.set_ylabel("fraction of instances")
    axes.set_title(title)
    axes.grid(True, which="both", alpha=0.3)
    axes.legend(loc="lower right")
    figure.tight_layout()
    figure.savefig(path)
    plt.close(figure)


def build_comparison(title, before, before_meta, after, after_meta, profile_path, changed_only, max_length=None,
                     embed_profile=None):
    """Before/after report; the profile is drawn to profile_path (if any) and embedded when embed_profile
    (default: whether it was drawn)."""
    all_names = sorted(set(before) | set(after), key=lambda name: (name.casefold(), name))
    common = [name for name in all_names if name in before and name in after]
    both_solved = [name for name in common if before[name].solved and after[name].solved]
    newly_solved = [name for name in common if after[name].solved and not before[name].solved]
    newly_failed = [name for name in common if before[name].solved and not after[name].solved]
    solved_before = sum(1 for name in common if before[name].solved)
    solved_after = sum(1 for name in common if after[name].solved)

    counted = [name for name in both_solved if has_counts(before[name]) and has_counts(after[name])]
    sgm_evaluations_before = shifted_geometric_mean([before[name].objective_evaluations for name in counted],
                                                    EVALUATIONS_SHIFT)
    sgm_evaluations_after = shifted_geometric_mean([after[name].objective_evaluations for name in counted],
                                                   EVALUATIONS_SHIFT)
    timed = [name for name in both_solved if before[name].solve is not None and after[name].solve is not None]
    sgm_time_before = shifted_geometric_mean([before[name].solve for name in timed], TIME_SHIFT)
    sgm_time_after = shifted_geometric_mean([after[name].solve for name in timed], TIME_SHIFT)

    label_before = f"before: PR #{before_meta.get('pr', '?')}"
    if before_meta.get("sha"):
        label_before += f" ({before_meta['sha'][:7]})"
    label_after = f"after: PR #{after_meta.get('pr', '?')}"
    if after_meta.get("sha"):
        label_after += f" ({after_meta['sha'][:7]})"
    profile_drawn = False
    if profile_path and common:
        taus, rho_before, rho_after = performance_profile(before, after, common)
        plot_profile(profile_path, taus, rho_before, rho_after, label_before, label_after,
                     f"{title}: objective evaluations ({len(common)} instances)")
        profile_drawn = True
    if embed_profile is None:
        embed_profile = profile_drawn

    lines = [
        f"## {title} — before/after",
        "",
        f"**Before:** {run_description(before_meta)}  ",
        f"**After:** {run_description(after_meta)} · {configuration(after_meta)}",
        "",
    ]
    for key, label in (("preset", "preset"), ("time_limit", "time limit"), ("selection", "instance selection"),
                       ("instances_tree", "CUTE instance archives")):
        if before_meta.get(key) is not None and after_meta.get(key) is not None \
                and before_meta[key] != after_meta[key]:
            lines += [f"> ⚠️ Different {label}: before `{before_meta[key]}`, after `{after_meta[key]}`.", ""]
    if before_meta.get("linear_solver") != after_meta.get("linear_solver"):
        lines += [f"> Linear solver: before `{before_meta.get('linear_solver')}`, "
                  f"after `{after_meta.get('linear_solver')}`.", ""]
    if embed_profile:
        lines += [f"![Performance profile (objective evaluations)]({PROFILE_URL_PLACEHOLDER})", ""]
    lines += [
        "| | before | after | after/before |",
        "|:--|--:|--:|--:|",
        f"| solved ({len(common)} common instances) | {solved_before} | "
        f"{solved_after} ({solved_after - solved_before:+d}) | — |",
        f"| shifted geo. mean objective evaluations, solved by both ({len(counted)}) | "
        f"{format_number(sgm_evaluations_before)} | {format_number(sgm_evaluations_after)} | "
        f"{format_ratio(sgm_evaluations_before, sgm_evaluations_after)} |",
        f"| shifted geo. mean solve time, solved by both ({len(timed)}) | {format_time(sgm_time_before)} | "
        f"{format_time(sgm_time_after)} | {format_ratio(sgm_time_before, sgm_time_after)} |",
        "",
        f"<sub>Performance profile: ratio of the objective evaluations to the fewest of both runs, "
        f"a failure counting as ∞. Shifts: {EVALUATIONS_SHIFT:g} evaluations, {TIME_SHIFT:g}s.</sub>",
        "",
    ]
    if newly_solved:
        lines += [f"🟢 **Newly solved ({len(newly_solved)}):** " + format_names(newly_solved), ""]
    if newly_failed:
        lines += [f"🔴 **Newly failed ({len(newly_failed)}):** " + format_names(newly_failed), ""]
    if not newly_solved and not newly_failed:
        lines += ["No change in solved/failed status.", ""]
    only_before = [name for name in all_names if name not in after]
    only_after = [name for name in all_names if name not in before]
    if only_before:
        lines += [f"Only in before run ({len(only_before)}): " + format_names(only_before), ""]
    if only_after:
        lines += [f"Only in after run ({len(only_after)}): " + format_names(only_after), ""]

    def cell(result):
        if result is None:
            return "—"
        icon = "✅" if result.solved else "❌"
        count = result.objective_evaluations
        return f"{icon} {count}" if count is not None else f"{icon} —"

    rows = []
    for name in all_names:
        b, a = before.get(name), after.get(name)
        flag, ratio = "", "—"
        if b is not None and a is not None:
            if a.solved and not b.solved:
                flag = "🟢"
            elif b.solved and not a.solved:
                flag = "🔴"
            elif has_counts(a) and has_counts(b):
                count_before, count_after = max(b.objective_evaluations, 1), max(a.objective_evaluations, 1)
                ratio = f"{count_after / count_before:.2f}×"
                if count_after / count_before >= EVALUATIONS_RATIO_THRESHOLD:
                    flag = "📈"
                elif count_before / count_after >= EVALUATIONS_RATIO_THRESHOLD:
                    flag = "📉"
        if changed_only and not flag:
            continue
        rows.append(f"| {flag} | `{name}` | {cell(b)} | {cell(a)} | {ratio} |")

    console = [f"Solved before: {solved_before} / {len(common)}, after: {solved_after} / {len(common)}",
               f"Shifted geometric mean of the objective evaluations ({len(counted)} solved by both): "
               f"before {format_number(sgm_evaluations_before)}, after {format_number(sgm_evaluations_after)}",
               f"Newly solved: {len(newly_solved)}, newly failed: {len(newly_failed)}"]
    if not rows:
        return "\n".join(lines), console, profile_drawn

    summary = (f"Instances whose status or objective evaluations changed ({len(rows)} of {len(all_names)})"
               if changed_only else f"Per-instance objective evaluations ({len(all_names)} instances)")
    header = [f"<details><summary>{summary}</summary>", "",
              "| | instance | before | after | after/before |", "|:-:|:--|--:|--:|--:|"]
    footer = ["", f"<sub>🟢 newly solved · 🔴 newly failed · 📉/📈 at least {EVALUATIONS_RATIO_THRESHOLD:g}× "
                  f"fewer/more objective evaluations (both solved)</sub>", "", "</details>", ""]
    report = "\n".join(lines + header + rows + footer)
    if max_length is not None and len(report) > max_length:
        # keep the report postable: drop table rows until it fits
        kept, note = list(rows), ""
        while kept and len("\n".join(lines + header + kept + [note] + footer)) > max_length:
            kept = kept[:max(0, len(kept) - max(1, len(kept) // 10))]
            note = f"| | … {len(rows) - len(kept)} more rows, see the run summary or the results artifact | | | |"
        report = "\n".join(lines + header + kept + [note] + footer)
    return report, console, profile_drawn


def set_output(name, value):
    path = os.environ.get("GITHUB_OUTPUT")
    if path:
        with open(path, "a") as file:
            file.write(f"{name}={value}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--title", default="CUTE benchmark")
    parser.add_argument("--results-dir", default="results")
    parser.add_argument("--output", default="benchmark-comment.md", help="report posted as the PR comment")
    parser.add_argument("--full-output", default=None,
                        help="optional: report with every instance, e.g. for the run summary")
    parser.add_argument("--compare-pr", default="", help="PR whose latest results are the 'before' (optional)")
    parser.add_argument("--compare-artifact", default="",
                        help="name of the results artifact of --compare-pr, e.g. cute-results-ipopt-pr-123")
    parser.add_argument("--profile", default="performance-profile.png", help="path of the performance profile")
    args = parser.parse_args()

    with open(os.path.join(args.results_dir, "results.tsv")) as file:
        after = parse_results_tsv(file.read())
    with open(os.path.join(args.results_dir, "meta.json")) as file:
        after_meta = json.load(file)
    if not after:
        print("::error::No results to report")
        return 1

    summary, console = build_summary(args.title, after, after_meta)

    fetched = None
    if args.compare_pr:
        try:
            fetched = fetch_from_artifact(os.environ["GITHUB_REPOSITORY"], args.compare_artifact,
                                          os.environ["GITHUB_TOKEN"], os.environ.get("GITHUB_RUN_ID", ""))
        except Exception as error:  # noqa: BLE001 - any failure falls back to the plain report
            print(f"Could not read the artifact {args.compare_artifact} ({error})")
        if fetched is None:
            print(f"::warning title={annotation_title(args.title)}::No results artifact {args.compare_artifact} "
                  f"found for PR #{args.compare_pr}; posting the plain report instead")

    if fetched is None:
        report, compared, drawn = summary, False, False
        full_report, _ = build_summary(args.title, after, after_meta, with_table=True)
    else:
        before, before_meta = fetched
        before_meta.setdefault("pr", args.compare_pr)
        print(f"Comparing against PR #{args.compare_pr} ({len(before)} instances)")
        report, console, drawn = build_comparison(args.title, before, before_meta, after, after_meta, args.profile,
                                                  changed_only=True, max_length=MAX_COMMENT_LENGTH)
        # the full report (run summary) embeds the same image
        full_report, _, _ = build_comparison(args.title, before, before_meta, after, after_meta, None,
                                             changed_only=False, embed_profile=drawn)
        compared = True

    with open(args.output, "w") as file:
        file.write(report)
    if args.full_output:
        with open(args.full_output, "w") as file:
            file.write(full_report)
    print("")
    print("================= SUMMARY =================")
    for line in console:
        print(line)
    print("===========================================")
    set_output("comparison", "true" if compared else "false")
    set_output("profile", "true" if drawn else "false")
    return 0


if __name__ == "__main__":
    sys.exit(main())
