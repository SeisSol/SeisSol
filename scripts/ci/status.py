#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
"""The state of the CI of a commit, a branch or a pull request, and why its jobs have failed.

Without logging in to GitHub, the logs of the jobs cannot be read; what can be read are the checks
of a commit, the steps of their jobs (which of them failed), and the artifacts of the runs, through
nightly.link. Failed jobs upload what tells why (see the workflows: build reports, the results of
the tests as JUnit with their output, the findings of clang-tidy); this script reads them and
prints the failed tests and findings, each with an excerpt of its output.

    scripts/ci/status.py davschneller/simcopy     # a branch
    scripts/ci/status.py 1615                     # a pull request
    scripts/ci/status.py 5eae63044                # a commit

GitHub answers 60 requests an hour without a token; a token in GITHUB_TOKEN or GH_TOKEN raises that,
fetches the artifacts from GitHub itself, and lets --logs print the end of the log of each failed
step. Only the standard library of Python is needed.
"""

import argparse
import io
import json
import os
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

API = "https://api.github.com"
USER_AGENT = "seissol-ci-status"

#: the artifacts that tell why a job failed, by the start of their names; the others (e.g. the
#: installed SeisSol, the coverage) are not downloaded
DIAGNOSTIC_ARTIFACTS = (
    "build-report-",
    "logs-",
    "sanitizer-",
    "aarch64-",
    "emulated-",
    "tidy-",
    "pytest-",
    "coverage-logs",
    "macos-logs",
    "windows-logs",
    "report-",
)

#: lines of the output of a test that tell what went wrong
ERROR_LINE = re.compile(
    r"error|fail|assert|abort|sanitizer|segmentation|signal|exception|timeout|deviat",
    re.IGNORECASE,
)
#: lines that say so without telling what: the verdict of a test, and the information of the log of
#: SeisSol (e.g. "Check the mesh for geometric errors"), unless a sanitizer writes into it
GENERIC_LINE = re.compile(
    r"^(Failed|error: exit code -?\d+)$|\binfo \d+ :(?!.*runtime error)"
)


class GitHubError(Exception):
    pass


class GitHub:
    """The requests to GitHub, with a token if there is one."""

    def __init__(self, repo, token=None):
        self.repo = repo
        self.token = token
        self.requests = 0
        self.remaining = None

    def _open(self, url, accept="application/vnd.github+json"):
        headers = {"User-Agent": USER_AGENT, "Accept": accept}
        if self.token is not None and url.startswith(API):
            headers["Authorization"] = f"Bearer {self.token}"
        request = urllib.request.Request(url, headers=headers)
        try:
            response = urllib.request.urlopen(request, timeout=120)
        except urllib.error.HTTPError as error:
            if (
                error.code in (403, 429)
                and error.headers.get("X-RateLimit-Remaining") == "0"
            ):
                raise GitHubError(
                    "GitHub's limit of requests without a token is reached (60 an hour); "
                    "set GITHUB_TOKEN, or wait."
                ) from error
            raise GitHubError(f"{error.code} {error.reason}: {url}") from error
        if url.startswith(API):
            self.requests += 1
            self.remaining = response.headers.get(
                "X-RateLimit-Remaining", self.remaining
            )
        return response

    def get(self, path):
        with self._open(f"{API}/repos/{self.repo}/{path}") as response:
            return json.load(response)

    def get_all(self, path, key):
        """All pages of a list in the field `key` of the answers."""
        items = []
        page = 1
        while True:
            separator = "&" if "?" in path else "?"
            answer = self.get(f"{path}{separator}per_page=100&page={page}")
            items += answer[key]
            if len(answer[key]) < 100:
                return items
            page += 1

    def artifact(self, artifact_id):
        """The zip file of an artifact: through nightly.link without a token."""
        if self.token is None:
            url = (
                f"https://nightly.link/{self.repo}/actions/artifacts/{artifact_id}.zip"
            )
        else:
            url = f"{API}/repos/{self.repo}/actions/artifacts/{artifact_id}/zip"
        with self._open(url, accept="*/*") as response:
            return response.read()

    def job_log(self, job_id):
        with self._open(
            f"{API}/repos/{self.repo}/actions/jobs/{job_id}/logs", accept="*/*"
        ) as response:
            return response.read().decode("utf-8", errors="replace")


def resolve(github, target):
    """The commit of a pull request (a number), a branch or a commit, and what it is."""
    if re.fullmatch(r"\d{1,6}", target):
        pull = github.get(f"pulls/{target}")
        return pull["head"]["sha"], f"pull request #{target} ({pull['head']['ref']})"
    commit = github.get(f"commits/{urllib.parse.quote(target, safe='')}")
    if commit["sha"].startswith(target):
        return commit["sha"], "commit"
    return commit["sha"], f"branch {target}"


def artifact_matches(job_name, artifact_name):
    """Whether an artifact belongs to a job, by the names the workflows give them."""
    tidy = re.fullmatch(r"clang-tidy(-gpu)? \((\d+)/\d+\)", job_name)
    if tidy:
        return artifact_name == f"tidy-{'gpu' if tidy[1] else 'cpu'}-{tidy[2]}"
    run = re.fullmatch(r"seissol-run-test<([^,>]+),.*><r(\d+)×t(\d+)>", job_name)
    if run:
        return artifact_name == f"report-{run[1]}-r{run[2]}t{run[3]}"
    if job_name.startswith("pytest-"):
        return artifact_name == job_name
    single = {
        "coverage": "coverage-logs",
        "build-macos": "macos-logs",
        "build-windows": "windows-logs",
    }
    if job_name in single:
        return artifact_name == single[job_name]
    # e.g. seissol-build-test<gcc-13, hsw, Release>, sanitizer<thread>, aarch64<neon>: an artifact
    # whose name has all parts of the matrix
    matrix = re.fullmatch(r"[\w-]+<(.+)>", job_name)
    if matrix:
        parts = [part.strip() for part in matrix[1].split(",")]
        return artifact_name.startswith(DIAGNOSTIC_ARTIFACTS) and all(
            part in artifact_name for part in parts
        )
    return False


def excerpt(text, lines):
    """The first lines of an output that tell what went wrong, or else its end; at most `lines`."""
    if not text:
        return []
    all_lines = [line.rstrip() for line in text.splitlines() if line.strip()]
    errors = [
        line
        for line in all_lines
        if ERROR_LINE.search(line) and not GENERIC_LINE.search(line.strip())
    ]
    return errors[:lines] if errors else all_lines[-lines:]


def junit_failures(content, lines):
    """The failed tests of a JUnit file (of ctest or pytest), with an excerpt of their output."""
    failures = []
    root = ET.fromstring(content)
    for case in root.iter("testcase"):
        problem = case.find("failure")
        if problem is None:
            problem = case.find("error")
        if problem is None and case.get("status") not in ("fail", "failed"):
            continue
        name = case.get("name", "?")
        if case.get("classname") and case.get("classname") not in name:
            name = f"{case.get('classname')}::{name}"
        text = ""
        if problem is not None:
            text = "\n".join(filter(None, [problem.get("message"), problem.text]))
        output = "\n".join(node.text or "" for node in case.iter("system-out"))
        failures.append((name, excerpt(f"{text}\n{output}", lines)))
    return failures


def build_report_failures(content):
    """What failed in the variants of a report of compile_variants.py."""
    report = json.loads(content)
    found = []
    for variant in report.get("variants", []):
        if variant.get("status") != "failed":
            continue
        testing = variant.get("testing") or {}
        tests = testing.get("failed_tests") or []
        detail = f"{variant['id']} failed in {variant.get('phase_failed')}"
        if testing.get("total"):
            detail += (
                f": {testing.get('failed', len(tests))} of {testing['total']} tests"
            )
        found.append((detail, tests))
    return found


def build_log_errors(content, lines):
    """The errors of the compilers and of CMake in a log of a build."""
    pattern = re.compile(
        r"error:|error [A-Z]+\d+|CMake Error|FAILED:|ninja: build stopped"
    )
    errors = [line.rstrip() for line in content.splitlines() if pattern.search(line)]
    return errors[:lines]


def read_artifact(name, data, lines):
    """What an artifact tells about a failure, as lines to print."""
    out = []
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        files = archive.namelist()
        junit_seen = False
        for file in sorted(files):
            base = file.rsplit("/", 1)[-1]
            content = archive.read(file)
            if base == "build-report.json":
                for detail, tests in build_report_failures(content):
                    out.append(f"{detail}")
                    out += [f"  - {test}" for test in tests]
            elif base.endswith(".xml") and not base.startswith("coverage"):
                try:
                    failures = junit_failures(content, lines)
                except ET.ParseError:
                    continue
                junit_seen = True
                if failures:
                    out.append(f"{file}: {len(failures)} failed")
                for test, text in failures:
                    out.append(f"  - {test}")
                    out += [f"      {line}" for line in text]
            elif base == "tidy-findings.txt":
                findings = content.decode("utf-8", errors="replace").splitlines()
                out.append(f"{len(findings)} findings of clang-tidy")
                out += [f"  {finding}" for finding in findings[: 5 * lines]]
            elif base == "build.log":
                errors = build_log_errors(
                    content.decode("utf-8", errors="replace"), lines
                )
                if errors:
                    out.append(f"{file}: errors of the build")
                    out += [f"  {line}" for line in errors]
            elif base == "report.txt":
                text = content.decode("utf-8", errors="replace")
                out.append(f"{file}:")
                out += [f"  {line}" for line in excerpt(text, lines)]
        if not junit_seen:
            for file in files:
                if file.endswith("LastTestsFailed.log"):
                    failed = (
                        archive.read(file).decode("utf-8", errors="replace").split()
                    )
                    out.append(f"{file}: {len(failed)} failed")
                    out += [f"  - {test}" for test in failed]
    return out


def run_of(check_run):
    """The workflow run of the check run of a job (whose id the check run shares)."""
    return int(re.search(r"/actions/runs/(\d+)/", check_run["html_url"])[1])


def failed_steps(job):
    return [
        step["name"]
        for step in job.get("steps", [])
        if step.get("conclusion") == "failure"
    ]


def log_excerpt(log, lines):
    """The end of the log of a job up to its first error."""
    all_lines = log.splitlines()
    for index, line in enumerate(all_lines):
        if "##[error]" in line:
            return all_lines[max(0, index - lines) : index + 1]
    return all_lines[-lines:]


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "target", help="a branch, a commit, or the number of a pull request"
    )
    parser.add_argument(
        "--repo", default="SeisSol/SeisSol", help="the repository (owner/name)"
    )
    parser.add_argument(
        "--lines",
        type=int,
        default=8,
        help="the lines of output per failure (default: 8)",
    )
    parser.add_argument(
        "--no-artifacts",
        action="store_true",
        help="only the checks and their failed steps",
    )
    parser.add_argument(
        "--logs",
        action="store_true",
        help="print the end of the log of each failed job (needs a token)",
    )
    args = parser.parse_args(argv)

    token = os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")
    github = GitHub(args.repo, token)
    try:
        sha, what = resolve(github, args.target)
        runs = github.get_all(f"commits/{sha}/check-runs", "check_runs")
        print(f"{args.repo} at {sha[:12]} ({what})")

        states = {}
        for run in runs:
            state = run["conclusion"] or run["status"]
            states[state] = states.get(state, 0) + 1
        print(
            "checks: "
            + ", ".join(f"{count} {state}" for state, count in sorted(states.items()))
        )

        failed = [run for run in runs if run["conclusion"] in ("failure", "timed_out")]
        pending = [run["name"] for run in runs if run["status"] != "completed"]
        if pending:
            print(f"still running or queued: {', '.join(sorted(pending))}")
        if not failed:
            print("no check has failed.")
            return 0

        # the jobs and the artifacts of the runs of the failed checks, a request or two per run
        jobs = {}
        artifacts = {}
        for run in failed:
            run_id = run_of(run)
            if run_id not in jobs:
                jobs[run_id] = {
                    job["id"]: job
                    for job in github.get_all(f"actions/runs/{run_id}/jobs", "jobs")
                }
                if not args.no_artifacts:
                    artifacts[run_id] = github.get_all(
                        f"actions/runs/{run_id}/artifacts", "artifacts"
                    )
        for run in sorted(failed, key=lambda run: run["name"]):
            print()
            run_id = run_of(run)
            steps = failed_steps(jobs[run_id].get(run["id"], {}))
            print(
                f"FAILED {run['name']}" + (f" (in {', '.join(steps)})" if steps else "")
            )
            print(f"  {run['html_url']}")
            if args.logs:
                if token is None:
                    print("  (the log needs a token in GITHUB_TOKEN)")
                else:
                    for line in log_excerpt(github.job_log(run["id"]), 3 * args.lines):
                        print(f"  | {line}")
            if args.no_artifacts:
                continue
            mine = [
                artifact
                for artifact in artifacts[run_id]
                if artifact_matches(run["name"], artifact["name"])
                and not artifact["expired"]
            ]
            if not mine:
                print("  (no artifact of it tells more)")
            for artifact in mine:
                print(f"  {artifact['name']}:")
                try:
                    found = read_artifact(
                        artifact["name"], github.artifact(artifact["id"]), args.lines
                    )
                except (GitHubError, zipfile.BadZipFile) as error:
                    print(f"    (could not read it: {error})")
                    continue
                for line in found or ["(nothing about a failure in it)"]:
                    print(f"    {line}")
    except GitHubError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    finally:
        if github.remaining is not None and int(github.remaining) < 15:
            print(
                f"\n({github.remaining} requests to GitHub left this hour)",
                file=sys.stderr,
            )
    return 1


if __name__ == "__main__":
    sys.exit(main())
