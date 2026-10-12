#!/usr/bin/env python3
# Copyright (c) 2026 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

"""Ctrl+C end-to-end driver (POSIX only).

Starts a child solve, waits until Uno has completed a few outer iterations, sends a
real SIGINT and checks that the child terminates normally and prints the expected
markers. The child prints, after the solve:
   STATUS <optimization status name>
   ITERATIONS <number of outer iterations>
plus any extra markers passed with --expect.

usage: driver.py [--min-iterations N] [--expect MARKER]... -- command [args...]
"""

import argparse
import re
import signal
import subprocess
import sys
import threading
import time

# rows of Uno's iteration table, e.g. " 12      OPT    5.0000e+01 ..."
ITERATION_ROW = re.compile(r"^\s*\d+\s+[A-Z]+\s")


def fail(message, lines):
    print("------ child output ------")
    print("".join(lines), end="")
    print(f"::error::{message}")
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-iterations", type=int, default=3, help="outer iterations to complete before Ctrl+C")
    parser.add_argument("--startup-timeout", type=float, default=600., help="seconds allowed to reach --min-iterations")
    parser.add_argument("--shutdown-timeout", type=float, default=60., help="seconds allowed to terminate after Ctrl+C")
    parser.add_argument("--expect", action="append", default=[], help="marker line the child must print")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    arguments = parser.parse_args()
    command = arguments.command[1:] if arguments.command[:1] == ["--"] else arguments.command
    if not command:
        parser.error("missing child command")

    # The child must start with the default SIGINT disposition, as in a terminal. A process
    # started in the background by a non-interactive shell inherits SIG_IGN instead, and then
    # e.g. Python does not install its KeyboardInterrupt handler.
    child = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1,
        preexec_fn=lambda: signal.signal(signal.SIGINT, signal.SIG_DFL))

    lines = []
    iteration_rows = 0
    lock = threading.Lock()

    def read_output():
        nonlocal iteration_rows
        for line in child.stdout:
            with lock:
                lines.append(line)
                if ITERATION_ROW.match(line):
                    iteration_rows += 1

    reader = threading.Thread(target=read_output, daemon=True)
    reader.start()

    def completed_iterations():
        with lock:
            return max(iteration_rows - 1, 0) # minus the initial point (iteration 0)

    # wait for a few outer iterations
    deadline = time.monotonic() + arguments.startup_timeout
    while child.poll() is None and completed_iterations() < arguments.min_iterations:
        if time.monotonic() > deadline:
            child.kill()
            fail(f"the child did not reach {arguments.min_iterations} iterations within {arguments.startup_timeout}s", lines)
        time.sleep(0.05)
    if child.poll() is not None:
        reader.join()
        fail(f"the child exited (code {child.returncode}) before Ctrl+C was sent: solved too quickly?", lines)

    iterations_at_signal = completed_iterations()
    print(f"Sending SIGINT after {iterations_at_signal} iterations", flush=True)
    child.send_signal(signal.SIGINT)

    try:
        child.wait(timeout=arguments.shutdown_timeout)
    except subprocess.TimeoutExpired:
        child.kill()
        reader.join()
        fail(f"the child ignored SIGINT: still running {arguments.shutdown_timeout}s after Ctrl+C", lines)
    reader.join()
    print("".join(lines), end="")

    # -2 (Python) or 130 (shell convention): killed by SIGINT, i.e. the signal was not caught
    if child.returncode != 0:
        fail(f"the child exited with code {child.returncode} (-2 means killed by SIGINT)", lines)
    output = "".join(lines)
    if not re.search(r"^STATUS USER_INTERRUPT$", output, re.MULTILINE):
        fail("the child did not report STATUS USER_INTERRUPT", lines)
    match = re.search(r"^ITERATIONS (\d+)$", output, re.MULTILINE)
    if match is None:
        fail("the child did not report its number of iterations", lines)
    final_iterations = int(match.group(1))
    # the iteration in progress when the signal arrived may complete (+1), plus a printing race (+1)
    if final_iterations > iterations_at_signal + 2:
        fail(f"the child ran {final_iterations} iterations, expected at most {iterations_at_signal + 2}", lines)
    for marker in arguments.expect:
        if not re.search(rf"^{re.escape(marker)}$", output, re.MULTILINE):
            fail(f"the child did not print '{marker}'", lines)

    print(f"::notice title=SIGINT::interrupted after {final_iterations} iterations "
          f"(signal sent after {iterations_at_signal})")


if __name__ == "__main__":
    main()
