#!/usr/bin/env python3
"""Before/after comparison of two Mittelmann benchmark runs of uno_ampl.

The "after" results are those of the current run (results/results.tsv and
results/meta.json, written by the benchmark step). The "before" results are
those of the latest benchmark run on another PR, looked up in this order:
  1. the `mittelmann-results-pr-<N>` artifact uploaded by that run;
  2. otherwise, the last report posted in that PR's benchmark comment (this
     covers runs that predate the results artifact).

The comparison is written as Markdown to --output. If the "before" results
cannot be found, a warning is printed, --output is left untouched and the
script exits with 0 so that the plain report is posted instead.

Standard library only.
"""

import argparse
import io
import json
import math
import os
import re
import sys
import urllib.error
import urllib.request
import zipfile

API_URL = os.environ.get("GITHUB_API_URL", "https://api.github.com")
ARTIFACT_PREFIX = "mittelmann-results-pr-"
# shift (in seconds) of the shifted geometric mean, as in Mittelmann's benchmarks
SGM_SHIFT = 10.0
# relative change in time (after/before) beyond which an instance is flagged
TIME_RATIO_THRESHOLD = 1.25
# instances faster than this (in both runs) are never flagged: timing noise
MIN_FLAGGED_TIME = 1.0


class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, *args, **kwargs):
        return None


def github_request(path, token):
    request = urllib.request.Request(API_URL + path, headers={
        "Authorization": f"Bearer {token}",
        "Accept": "application/vnd.github+json",
        "X-GitHub-Api-Version": "2022-11-28",
    })
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.load(response)


def download_artifact_zip(repository, artifact_id, token):
    # the API answers with a redirect to a signed URL; the token must not be sent there
    request = urllib.request.Request(f"{API_URL}/repos/{repository}/actions/artifacts/{artifact_id}/zip", headers={
        "Authorization": f"Bearer {token}",
        "Accept": "application/vnd.github+json",
    })
    opener = urllib.request.build_opener(NoRedirect)
    try:
        with opener.open(request, timeout=60) as response:
            return response.read()
    except urllib.error.HTTPError as error:
        if error.code not in (301, 302, 303, 307, 308):
            raise
        location = error.headers["Location"]
    with urllib.request.urlopen(location, timeout=120) as response:
        return response.read()


def parse_results_tsv(text):
    results = {}
    for line in text.splitlines():
        fields = line.strip().split("\t")
        if len(fields) != 3:
            continue
        name, status, seconds = fields
        results[name] = (status == "solved", float(seconds))
    return results


def fetch_from_artifact(repository, pr_number, token, current_run_id):
    name = f"{ARTIFACT_PREFIX}{pr_number}"
    listing = github_request(f"/repos/{repository}/actions/artifacts?name={name}&per_page=100", token)
    artifacts = [artifact for artifact in listing.get("artifacts", [])
                 if not artifact.get("expired")
                 and str((artifact.get("workflow_run") or {}).get("id")) != str(current_run_id)]
    if not artifacts:
        return None
    latest = max(artifacts, key=lambda artifact: artifact["created_at"])
    with zipfile.ZipFile(io.BytesIO(download_artifact_zip(repository, latest["id"], token))) as archive:
        results = parse_results_tsv(archive.read("results.tsv").decode())
        meta = json.loads(archive.read("meta.json").decode()) if "meta.json" in archive.namelist() else {}
    if not results:
        return None
    meta.setdefault("source", "artifact")
    return results, meta


# rows of the plain report: | ✅ | `name` | solved in 12.34s |
ROW_PATTERN = re.compile(r"^\|\s*(✅|❌)\s*\|\s*`([^`]+)`\s*\|\s*(?:solved|failed) in ([0-9.]+)s\s*\|\s*$")
REPORT_HEADER = "## Mittelmann benchmark"


def fetch_from_comments(repository, pr_number, token):
    comments = []
    page = 1
    while True:
        batch = github_request(f"/repos/{repository}/issues/{pr_number}/comments?per_page=100&page={page}", token)
        comments.extend(batch)
        if len(batch) < 100:
            break
        page += 1
    # most recent report first: last comment, last appended section
    for comment in reversed(comments):
        body = comment.get("body") or ""
        if REPORT_HEADER not in body:
            continue
        for section in reversed(body.split(REPORT_HEADER)[1:]):
            results = {}
            for line in section.splitlines():
                match = ROW_PATTERN.match(line.strip())
                if match:
                    results[match.group(2)] = (match.group(1) == "✅", float(match.group(3)))
            if not results:
                continue
            meta = {"source": "comment", "comment_url": comment.get("html_url")}
            sha = re.search(r"Commit `([0-9a-f]{7,40})`", section)
            if sha:
                meta["sha"] = sha.group(1)
            run_url = re.search(r"\[run details\]\(([^)]+)\)", section)
            if run_url:
                meta["run_url"] = run_url.group(1)
            time_limit = re.search(r"time_limit=([0-9]+)s", section)
            if time_limit:
                meta["time_limit"] = int(time_limit.group(1))
            return results, meta
    return None


def shifted_geometric_mean(times, shift=SGM_SHIFT):
    if not times:
        return float("nan")
    return math.exp(sum(math.log(t + shift) for t in times) / len(times)) - shift


def describe_run(label, pr_number, meta):
    parts = [f"PR #{pr_number}"]
    if meta.get("sha"):
        parts.append(f"commit `{meta['sha'][:7]}`")
    if meta.get("run_url"):
        parts.append(f"[run]({meta['run_url']})")
    elif meta.get("comment_url"):
        parts.append(f"[report]({meta['comment_url']})")
    return f"**{label}:** " + " · ".join(parts)


def format_cell(result):
    if result is None:
        return "—"
    solved, seconds = result
    return f"{'✅' if solved else '❌'} {seconds:.2f}s"


def build_report(before, before_meta, before_pr, after, after_meta, after_pr):
    time_limit = after_meta.get("time_limit")
    all_names = sorted(set(before) | set(after))
    common = [name for name in all_names if name in before and name in after]

    newly_solved = [name for name in common if after[name][0] and not before[name][0]]
    newly_failed = [name for name in common if before[name][0] and not after[name][0]]
    solved_before = sum(1 for name in common if before[name][0])
    solved_after = sum(1 for name in common if after[name][0])

    # shifted geometric means over the common instances, failures counted at the time limit
    def capped(result):
        solved, seconds = result
        return seconds if solved or time_limit is None else max(seconds, float(time_limit))
    sgm_before = shifted_geometric_mean([capped(before[name]) for name in common])
    sgm_after = shifted_geometric_mean([capped(after[name]) for name in common])
    both_solved = [name for name in common if before[name][0] and after[name][0]]
    sgm_both_before = shifted_geometric_mean([before[name][1] for name in both_solved])
    sgm_both_after = shifted_geometric_mean([after[name][1] for name in both_solved])

    def ratio(a, b):
        return f"{b / a:.2f}×" if a > 0 and not math.isnan(a) and not math.isnan(b) else "—"

    title = "## Mittelmann benchmark — before/after"
    if time_limit is not None:
        title += f" (time_limit={time_limit}s)"
    lines = [
        title,
        "",
        describe_run("Before", before_pr, before_meta) + "  ",
        describe_run("After", after_pr, after_meta),
        "",
    ]
    if before_meta.get("time_limit") not in (None, time_limit):
        lines += [f"> ⚠️ The before run used time_limit={before_meta['time_limit']}s.", ""]
    lines += [
        f"| | before | after | after/before |",
        f"|:--|--:|--:|--:|",
        f"| solved ({len(common)} common instances) | {solved_before} | {solved_after} | {solved_after - solved_before:+d} |",
        f"| shifted geo. mean time, all common instances¹ | {sgm_before:.2f}s | {sgm_after:.2f}s | {ratio(sgm_before, sgm_after)} |",
        f"| shifted geo. mean time, solved by both ({len(both_solved)}) | {sgm_both_before:.2f}s | {sgm_both_after:.2f}s | {ratio(sgm_both_before, sgm_both_after)} |",
        "",
        f"<sub>¹ shift of {SGM_SHIFT:g}s; failures counted at the time limit.</sub>",
        "",
    ]
    if newly_solved:
        lines.append(f"🟢 **Newly solved ({len(newly_solved)}):** " + ", ".join(f"`{name}`" for name in newly_solved))
        lines.append("")
    if newly_failed:
        lines.append(f"🔴 **Newly failed ({len(newly_failed)}):** " + ", ".join(f"`{name}`" for name in newly_failed))
        lines.append("")
    if not newly_solved and not newly_failed:
        lines += ["No change in solved/failed status.", ""]
    only_before = [name for name in all_names if name not in after]
    only_after = [name for name in all_names if name not in before]
    if only_before:
        lines += [f"Only in before run ({len(only_before)}): " + ", ".join(f"`{name}`" for name in only_before), ""]
    if only_after:
        lines += [f"Only in after run ({len(only_after)}): " + ", ".join(f"`{name}`" for name in only_after), ""]

    lines += [
        f"<details><summary>Per-instance results ({len(all_names)} instances)</summary>",
        "",
        "| | instance | before | after | after/before |",
        "|:-:|:--|--:|--:|--:|",
    ]
    for name in all_names:
        b, a = before.get(name), after.get(name)
        flag, time_ratio = "", "—"
        if b is not None and a is not None:
            if a[0] and not b[0]:
                flag = "🟢"
            elif b[0] and not a[0]:
                flag = "🔴"
            elif a[0] and b[0]:
                time_ratio = f"{a[1] / b[1]:.2f}×" if b[1] > 0 else "—"
                if max(a[1], b[1]) >= MIN_FLAGGED_TIME and b[1] > 0:
                    if a[1] / b[1] >= TIME_RATIO_THRESHOLD:
                        flag = "🐢"
                    elif b[1] / max(a[1], 1e-9) >= TIME_RATIO_THRESHOLD:
                        flag = "⚡"
        lines.append(f"| {flag} | `{name}` | {format_cell(b)} | {format_cell(a)} | {time_ratio} |")
    lines += [
        "",
        f"<sub>🟢 newly solved · 🔴 newly failed · ⚡/🐢 at least {TIME_RATIO_THRESHOLD:g}× faster/slower "
        f"(both solved, ≥ {MIN_FLAGGED_TIME:g}s)</sub>",
        "",
        "</details>",
        "",
    ]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compare-pr", required=True, help="PR whose latest results are the 'before'")
    parser.add_argument("--current-pr", required=True, help="PR of the current run (the 'after')")
    parser.add_argument("--results-dir", default="results")
    parser.add_argument("--output", default="benchmark-comment.md")
    args = parser.parse_args()

    repository = os.environ["GITHUB_REPOSITORY"]
    token = os.environ["GITHUB_TOKEN"]
    run_id = os.environ.get("GITHUB_RUN_ID", "")

    with open(os.path.join(args.results_dir, "results.tsv")) as file:
        after = parse_results_tsv(file.read())
    with open(os.path.join(args.results_dir, "meta.json")) as file:
        after_meta = json.load(file)

    fetched = None
    try:
        fetched = fetch_from_artifact(repository, args.compare_pr, token, run_id)
        if fetched is None:
            print(f"No results artifact found for PR #{args.compare_pr}; looking at its comments")
    except Exception as error:  # noqa: BLE001 - fall back to the comments
        print(f"Could not read the results artifact of PR #{args.compare_pr} ({error}); looking at its comments")
    if fetched is None:
        try:
            fetched = fetch_from_comments(repository, args.compare_pr, token)
        except Exception as error:  # noqa: BLE001 - any failure falls back to the plain report
            print(f"::warning title=Mittelmann comparison::Could not fetch results of PR #{args.compare_pr}: {error}")
            return 0
    if fetched is None:
        print(f"::warning title=Mittelmann comparison::No benchmark results found for PR #{args.compare_pr}; "
              "posting the plain report instead")
        return 0

    before, before_meta = fetched
    print(f"Comparing against PR #{args.compare_pr} ({len(before)} instances, from {before_meta['source']})")
    report = build_report(before, before_meta, args.compare_pr, after, after_meta, args.current_pr)
    with open(args.output, "w") as file:
        file.write(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
