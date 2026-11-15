# Copyright (c) 2026 Alexis Montoison and Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

# Add the libkrylov bundle (library + Julia runtime, see download_dependencies.sh) to repaired unopy wheels, under
# unopy/libkrylov. This runs after auditwheel/delocate/delvewheel so that the repair tools never rewrite the
# self-contained bundle (whose libraries locate each other via $ORIGIN/@loader_path and their own directory).
#
# Wheels cannot store symbolic links, and dereferencing all of them would duplicate most libraries. A library is
# only looked up under the names referenced in the binaries of the bundle (SONAMEs, install names, dlopen names), so
# each file is stored under its referenced names only (or under its own name if it is not referenced).
#
# usage: python add_libkrylov_to_wheel.py <libkrylov bundle directory> <wheel or directory of wheels>...

import base64
import hashlib
import os
import re
import sys
import tempfile
import zipfile

PACKAGE_DIRECTORY = "unopy/libkrylov"


def bundle_entries(bundle_directory):
    """Return the list of (relative path in the bundle, source file) to store in the wheel."""
    paths = []
    for root, _, files in os.walk(bundle_directory):
        for file in files:
            paths.append(os.path.join(root, file))
    # group the paths (files and symbolic links) by target file
    groups = {}
    for path in paths:
        groups.setdefault(os.path.realpath(path), []).append(path)
    # names referenced in the binaries: delimited by characters that cannot appear in a library name (NUL, path
    # separators, ':' in the preload list of libjulia, quotes, @rpath/...)
    contents = b"".join(open(target, "rb").read() for target in groups)
    def is_string_referenced(name):
        pattern = rb"(?<![A-Za-z0-9._+-])" + re.escape(name.encode()) + rb"(?![A-Za-z0-9._+-])"
        return re.search(pattern, contents) is not None

    def is_referenced(name):
        # Julia also opens libraries by stem + platform extension (e.g. "libopenlibm" -> libopenlibm.so)
        for extension in (".so", ".dylib", ".dll"):
            if name.endswith(extension) and is_string_referenced(name[:-len(extension)]):
                return True
        return is_string_referenced(name)

    entries = []
    for target, members in groups.items():
        kept = [path for path in members if is_referenced(os.path.basename(path))]
        if not kept:
            # the file itself (not one of its symbolic links)
            kept = [path for path in members if not os.path.islink(path)] or members[:1]
        for path in kept:
            entries.append((os.path.relpath(path, bundle_directory).replace(os.sep, "/"), target))
    return sorted(entries)


def record_hash(data):
    digest = base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b"=").decode()
    return f"sha256={digest}"


def add_bundle_to_wheel(wheel, entries):
    with zipfile.ZipFile(wheel) as source:
        names = source.namelist()
        record_name = next(name for name in names if name.endswith(".dist-info/RECORD"))
        record_lines = [line for line in source.read(record_name).decode().splitlines() if line]
        if any(name.startswith(PACKAGE_DIRECTORY + "/") for name in names):
            print(f"{wheel}: libkrylov is already bundled")
            return
        file_descriptor, temporary_wheel = tempfile.mkstemp(suffix=".whl", dir=os.path.dirname(os.path.abspath(wheel)))
        os.close(file_descriptor)
        with zipfile.ZipFile(temporary_wheel, "w", compression=zipfile.ZIP_DEFLATED) as destination:
            for info in source.infolist():
                if info.filename != record_name:
                    destination.writestr(info, source.read(info.filename))
            size = 0
            for relative_path, target in entries:
                arcname = f"{PACKAGE_DIRECTORY}/{relative_path}"
                with open(target, "rb") as file:
                    data = file.read()
                info = zipfile.ZipInfo(arcname, date_time=(1980, 1, 1, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = (os.stat(target).st_mode & 0o777 | 0o100000) << 16
                destination.writestr(info, data)
                record_lines.append(f"{arcname},{record_hash(data)},{len(data)}")
                size += len(data)
            record_lines = [line for line in record_lines if not line.startswith(record_name + ",")]
            record_lines.append(f"{record_name},,")
            destination.writestr(source.getinfo(record_name), "\n".join(record_lines) + "\n")
    os.chmod(temporary_wheel, os.stat(wheel).st_mode & 0o777)
    os.replace(temporary_wheel, wheel)
    print(f"{wheel}: added {len(entries)} libkrylov files ({size / 1e6:.0f} MB uncompressed)")


def main():
    if len(sys.argv) < 3:
        sys.exit("usage: python add_libkrylov_to_wheel.py <libkrylov bundle directory> <wheel or directory of wheels>...")
    bundle_directory = sys.argv[1]
    if not os.path.isdir(bundle_directory):
        # no libkrylov artifact for this platform (see download_dependencies.sh)
        print(f"{bundle_directory} does not exist: libkrylov is not bundled")
        return
    entries = bundle_entries(bundle_directory)
    wheels = []
    for argument in sys.argv[2:]:
        if os.path.isdir(argument):
            wheels += [os.path.join(argument, name) for name in os.listdir(argument) if name.endswith(".whl")]
        else:
            wheels.append(argument)
    if not wheels:
        sys.exit("no wheel found")
    for wheel in wheels:
        add_bundle_to_wheel(wheel, entries)


if __name__ == "__main__":
    main()
