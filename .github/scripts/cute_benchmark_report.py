#!/usr/bin/env python3
"""Report of a CUTE benchmark run of Uno, optionally compared against another PR.

The results of the current run are read from <results-dir>/results.tsv and
<results-dir>/meta.json (written by the CUTE benchmark workflow). Each row of
results.tsv is
    <name> <solved|failed> <wall seconds> <model read seconds> <solve seconds> <objective evaluations> <iterations>
the last four fields being empty when unknown (no optimization summary in the log).

Without --compare-pr, the report gives the number of solved instances and the
shifted geometric means of the objective evaluations and of the solve time and,
with a --reference, compares the run against the reference solver, with a
performance profile (reference-profile.png) with respect to the objective evaluations.

With --compare-pr, the "before" results are those of the latest
`<artifact-prefix><compare-pr>` artifact (another run of the same preset). The
report then compares both runs on their common instances, with Dolan-Moré
performance profiles with respect to the number of objective evaluations,
written to --profile-dir:
  - performance-profile.png: before vs after;
  - reference-profile.png: before, after and the --reference solver (if given),
    e.g. IPOPT for the ipopt preset, on the instances common to the three.
Matplotlib is needed as soon as a profile is drawn.
The report embeds each image as ![...](PLOT:<file name>); the workflow replaces
"PLOT:" with the URL of the published images. If the "before" results cannot be
found, a warning is printed and the plain report is written instead.

--reference is a table in the results.tsv format (times may be empty) whose
first line is "# Reference results ...: <solver label>".

"comparison=true|false" and "plots=true|false" (whether a profile was drawn) are
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

# images are embedded as ![...](PLOT:<file name>), "PLOT:" being replaced by the URL of the published images
PLOT_PREFIX = "PLOT:"
PROFILE_FILE = "performance-profile.png"
REFERENCE_PROFILE_FILE = "reference-profile.png"
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
    wall: Optional[float]
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
        if line.startswith("#") or len(fields) != 7:
            continue
        name, status, wall, read, solve, evaluations, iterations = fields
        results[name.strip()] = Result(status.strip() == "solved", optional(wall, float), optional(read, float),
                                       optional(solve, float), optional(evaluations, int), optional(iterations, int))
    return results


def read_reference(path):
    """(results, solver label) of a reference table."""
    with open(path) as file:
        text = file.read()
    first_line = text.splitlines()[0] if text else ""
    label = first_line.split(":", 1)[1].strip() if first_line.startswith("#") and ":" in first_line \
        else os.path.splitext(os.path.basename(path))[0]
    return parse_results_tsv(text), label


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


def build_summary(title, results, meta, with_table=False, reference=None, reference_label=None, images=()):
    """Report of a single run, against the reference solver if any (with_table: one row per instance,
    for the run summary)."""
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
    if reference is not None:
        lines += reference_section([results], ["this PR"], reference, reference_label, images)
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


def performance_profile(runs, names):
    """Dolan-Moré profile on the objective evaluations of several runs: (taus, [rho of each run]).
    A failed instance has an infinite ratio; an instance failed by all runs is never reached."""
    def cost(result):
        return max(result.objective_evaluations, 1) if has_counts(result) else math.inf

    ratios = [[] for _ in runs]
    for name in names:
        costs = [cost(run[name]) for run in runs]
        best = min(costs)
        for run_ratios, run_cost in zip(ratios, costs):
            run_ratios.append(run_cost / best if best < math.inf else math.inf)
    finite = [ratio for run_ratios in ratios for ratio in run_ratios if ratio < math.inf]
    tau_max = max(2.0, 1.1 * max(finite, default=1.0))
    taus = sorted(set([1.0, tau_max] + finite))
    return taus, [[sum(1 for ratio in run_ratios if ratio <= tau) / len(run_ratios) for tau in taus]
                  for run_ratios in ratios]


# one style per role, the same in every profile
BEFORE_STYLE = dict(color="#4269d0", linestyle="-")
AFTER_STYLE = dict(color="#ef8a0c", linestyle="--")  # also this PR, without comparison
REFERENCE_STYLE = dict(color="#3ca951", linestyle=":")


def plot_profile(path, taus, rhos, labels, styles, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    figure, axes = plt.subplots(figsize=(8, 5), dpi=120)
    for rho, label, style in zip(rhos, labels, styles):
        axes.step(taus, rho, where="post", label=label, linewidth=2, **style)
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


def run_label(prefix, meta):
    label = (f"{prefix}: " if prefix else "") + f"PR #{meta.get('pr', '?')}"
    return label + (f" ({meta['sha'][:7]})" if meta.get("sha") else "")


def common_names(*runs):
    return sorted(set.intersection(*(set(run) for run in runs)), key=lambda name: (name.casefold(), name))


def draw_profiles(directory, title, runs, labels, reference=None, reference_label=None):
    """Draws in directory the profile of the runs (if there are several: before/after) and, with a reference,
    the profile of the runs and the reference. Returns the file names of the drawn profiles."""
    drawn = []
    styles = [BEFORE_STYLE, AFTER_STYLE] if len(runs) > 1 else [AFTER_STYLE]
    names = common_names(*runs)
    if len(runs) > 1 and names:
        taus, rhos = performance_profile(runs, names)
        plot_profile(os.path.join(directory, PROFILE_FILE), taus, rhos, labels, styles,
                     f"{title}: objective evaluations ({len(names)} instances)")
        drawn.append(PROFILE_FILE)
    if reference is not None:
        names = common_names(*runs, reference)
        if names:
            taus, rhos = performance_profile([*runs, reference], names)
            plot_profile(os.path.join(directory, REFERENCE_PROFILE_FILE), taus, rhos, [*labels, reference_label],
                         [*styles, REFERENCE_STYLE], f"{title} vs {reference_label} ({len(names)} instances)")
            drawn.append(REFERENCE_PROFILE_FILE)
    return drawn


def reference_section(runs, headers, reference, reference_label, images):
    """The runs (this PR, or before and after) against the reference solver, on their common instances."""
    names = common_names(*runs, reference)
    if not names:
        return []
    columns = [*runs, reference]
    solved = [sum(1 for name in names if column[name].solved) for column in columns]
    counted = [name for name in names if all(has_counts(column[name]) for column in columns)]
    sgms = [shifted_geometric_mean([column[name].objective_evaluations for name in counted], EVALUATIONS_SHIFT)
            for column in columns]
    solved_by = "both" if len(columns) == 2 else f"the {len(columns)}"
    lines = [f"### Against {reference_label}", ""]
    if REFERENCE_PROFILE_FILE in images:
        lines += [f"![Performance profile against {reference_label} (objective evaluations)]"
                  f"({PLOT_PREFIX}{REFERENCE_PROFILE_FILE})", ""]
    lines += [
        "| | " + " | ".join([*headers, reference_label]) + " |",
        "|:--|" + "--:|" * len(columns),
        f"| solved ({len(names)} common instances) | " + " | ".join(str(count) for count in solved) + " |",
        f"| shifted geo. mean objective evaluations, solved by {solved_by} ({len(counted)}) | "
        + " | ".join(format_number(sgm) for sgm in sgms) + " |",
        "",
        f"<sub>{reference_label}: reference results from "
        f"[nonlinear_optimization_solver_benchmark](https://github.com/cvanaret/nonlinear_optimization_solver_benchmark), "
        f"obtained with its own tolerances and limits (times not comparable).</sub>",
        "",
    ]
    return lines


def build_comparison(title, before, before_meta, after, after_meta, images, changed_only, max_length=None,
                     reference=None, reference_label=None):
    """Before/after report (and against the reference solver, if any), embedding the profiles in images."""
    all_names = sorted(set(before) | set(after), key=lambda name: (name.casefold(), name))
    common = [name for name in all_names if name in before and name in after]
    both_solved = [name for name in common if before[name].solved and after[name].solved]
    newly_solved = [name for name in common if after[name].solved and not before[name].solved]
    newly_failed = [name for name in common if before[name].solved and not after[name].solved]
    solved_before = sum(1 for name in common if before[name].solved)
    solved_after = sum(1 for name in common if after[name].solved)

    counted = [name for name in both_solved if has_counts(before[name]) and has_counts(after[name])]
    different_counts = [name for name in counted
                        if before[name].objective_evaluations != after[name].objective_evaluations]
    fewer = sum(1 for name in different_counts if after[name].objective_evaluations < before[name].objective_evaluations)
    sgm_evaluations_before = shifted_geometric_mean([before[name].objective_evaluations for name in counted],
                                                    EVALUATIONS_SHIFT)
    sgm_evaluations_after = shifted_geometric_mean([after[name].objective_evaluations for name in counted],
                                                   EVALUATIONS_SHIFT)
    timed = [name for name in both_solved if before[name].solve is not None and after[name].solve is not None]
    sgm_time_before = shifted_geometric_mean([before[name].solve for name in timed], TIME_SHIFT)
    sgm_time_after = shifted_geometric_mean([after[name].solve for name in timed], TIME_SHIFT)

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
    if PROFILE_FILE in images:
        lines += [f"![Performance profile (objective evaluations)]({PLOT_PREFIX}{PROFILE_FILE})", ""]
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
    if different_counts:
        lines += [f"**Objective evaluations changed on {len(different_counts)} of the {len(counted)} instances "
                  f"solved by both:** fewer on {fewer}, more on {len(different_counts) - fewer}.", ""]
    else:
        lines += [f"Identical objective evaluations on all {len(counted)} instances solved by both.", ""]
    only_before = [name for name in all_names if name not in after]
    only_after = [name for name in all_names if name not in before]
    if only_before:
        lines += [f"Only in before run ({len(only_before)}): " + format_names(only_before), ""]
    if only_after:
        lines += [f"Only in after run ({len(only_after)}): " + format_names(only_after), ""]
    if reference is not None:
        lines += reference_section([before, after], ["before", "after"], reference, reference_label, images)

    def cell(result):
        if result is None:
            return "—"
        icon = "✅" if result.solved else "❌"
        count = result.objective_evaluations
        return f"{icon} {count}" if count is not None else f"{icon} —"

    rows = []
    for name in all_names:
        b, a = before.get(name), after.get(name)
        flag, ratio, changed = "", "—", False
        if b is not None and a is not None:
            if a.solved and not b.solved:
                flag, changed = "🟢", True
            elif b.solved and not a.solved:
                flag, changed = "🔴", True
            elif has_counts(a) and has_counts(b):
                count_before, count_after = max(b.objective_evaluations, 1), max(a.objective_evaluations, 1)
                ratio = f"{count_after / count_before:.2f}×"
                # any change in the number of objective evaluations is listed; large ones are flagged
                changed = b.objective_evaluations != a.objective_evaluations
                if count_after / count_before >= EVALUATIONS_RATIO_THRESHOLD:
                    flag = "📈"
                elif count_before / count_after >= EVALUATIONS_RATIO_THRESHOLD:
                    flag = "📉"
        if changed_only and not changed:
            continue
        rows.append(f"| {flag} | `{name}` | {cell(b)} | {cell(a)} | {ratio} |")

    console = [f"Solved before: {solved_before} / {len(common)}, after: {solved_after} / {len(common)}",
               f"Shifted geometric mean of the objective evaluations ({len(counted)} solved by both): "
               f"before {format_number(sgm_evaluations_before)}, after {format_number(sgm_evaluations_after)}",
               f"Newly solved: {len(newly_solved)}, newly failed: {len(newly_failed)}",
               f"Different objective evaluations (solved by both): {len(different_counts)} of {len(counted)}"]
    if not rows:
        return "\n".join(lines), console

    summary = (f"Instances whose status or objective evaluations changed ({len(rows)} of {len(all_names)})"
               if changed_only else f"Per-instance objective evaluations ({len(all_names)} instances)")
    header = [f"<details><summary>{summary}</summary>", "",
              "| | instance | before | after | after/before |", "|:-:|:--|--:|--:|--:|"]
    footer = ["", f"<sub>🟢 newly solved · 🔴 newly failed · 📉/📈 at least {EVALUATIONS_RATIO_THRESHOLD:g}× "
                  f"fewer/more objective evaluations (both solved); unflagged rows changed by less</sub>",
              "", "</details>", ""]
    report = "\n".join(lines + header + rows + footer)
    if max_length is not None and len(report) > max_length:
        # keep the report postable: drop table rows until it fits
        kept, note = list(rows), ""
        while kept and len("\n".join(lines + header + kept + [note] + footer)) > max_length:
            kept = kept[:max(0, len(kept) - max(1, len(kept) // 10))]
            note = f"| | … {len(rows) - len(kept)} more rows, see the run summary or the results artifact | | | |"
        report = "\n".join(lines + header + kept + [note] + footer)
    return report, console


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
    parser.add_argument("--profile-dir", default=".", help="directory of the performance profiles")
    parser.add_argument("--reference", default="", help="reference results of another solver (optional)")
    args = parser.parse_args()

    with open(os.path.join(args.results_dir, "results.tsv")) as file:
        after = parse_results_tsv(file.read())
    with open(os.path.join(args.results_dir, "meta.json")) as file:
        after_meta = json.load(file)
    if not after:
        print("::error::No results to report")
        return 1

    _, console = build_summary(args.title, after, after_meta)

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

    reference, reference_label = (None, None)
    if args.reference:
        reference, reference_label = read_reference(args.reference)
        print(f"Reference: {reference_label} ({len(reference)} instances)")

    if fetched is None:
        compared = False
        drawn = draw_profiles(args.profile_dir, args.title, [after], [run_label(None, after_meta)],
                              reference, reference_label)
        report, console = build_summary(args.title, after, after_meta, reference=reference,
                                        reference_label=reference_label, images=drawn)
        full_report, _ = build_summary(args.title, after, after_meta, with_table=True, reference=reference,
                                       reference_label=reference_label, images=drawn)
    else:
        before, before_meta = fetched
        before_meta.setdefault("pr", args.compare_pr)
        print(f"Comparing against PR #{args.compare_pr} ({len(before)} instances)")
        drawn = draw_profiles(args.profile_dir, args.title, [before, after],
                              [run_label("before", before_meta), run_label("after", after_meta)],
                              reference, reference_label)
        report, console = build_comparison(args.title, before, before_meta, after, after_meta, drawn,
                                           changed_only=True, max_length=MAX_COMMENT_LENGTH,
                                           reference=reference, reference_label=reference_label)
        full_report, _ = build_comparison(args.title, before, before_meta, after, after_meta, drawn,
                                          changed_only=False, reference=reference, reference_label=reference_label)
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
    set_output("plots", "true" if drawn else "false")
    return 0


if __name__ == "__main__":
    sys.exit(main())
