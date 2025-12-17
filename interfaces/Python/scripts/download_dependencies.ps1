Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

# Destination directory for third-party static libs
$dest = "interfaces/Python/dependencies/bqpd"
New-Item -ItemType Directory -Force -Path $dest | Out-Null
Set-Location $dest

# ---------------------------------

# --- Detect architecture ---
$arch = $env:PROCESSOR_ARCHITECTURE

# If running under WOW64 (x86 on x64), detect real architecture
if ($arch -eq "x86" -and $env:PROCESSOR_ARCHITEW6432) {
    $arch = $env:PROCESSOR_ARCHITEW6432
}

switch ($arch.ToLower()) {
    "amd64" { $ARCH = "x86_64" }
    "arm64" { $ARCH = "aarch64" }
    default {
        Write-Error "Unknown or unsupported architecture '$arch'. Modify this script."
        exit 1
    }
}

Write-Output "Detected architecture: $ARCH"

# --- Construct archive URL ---
$version = "v1.0.0"
$OS = "w64-mingw32"
$archiveName = "BQPD.$version.$ARCH-$OS-libgfortran5.tar.gz"
$assetUrl = "https://github.com/leyffer/BQPD_jll.jl/releases/download/BQPD-$version%2B0/$archiveName"
Write-Output "Downloading archive: $assetUrl"

Invoke-WebRequest -Uri $assetUrl -OutFile $archiveName

# ---------- Extract Archive ----------
Write-Output "Extracting tar.gz archive..."
tar -xzf $archiveName

# ---------- Show results ----------
Write-Output "Contents extracted:"
Get-ChildItem -Force

Write-Output "`nExpected result: libbqpd.a should now be present in $dest/lib"