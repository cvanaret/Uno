Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

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
$OS = "w64-mingw32"
Write-Output "Detected architecture: $ARCH $OS"

# destination directory for third-party static libs
$third_party_dir = "interfaces/Python/third_party"
New-Item -ItemType Directory -Force -Path $third_party_dir | Out-Null
Set-Location $third_party_dir

# download BQPD
$version = "v1.0.0"
$repo = "https://github.com/leyffer/BQPD_jll.jl/releases/download/BQPD-$version%2B0"
$asset_name = "BQPD.$version.$ARCH-$OS-libgfortran5.tar.gz"
$asset_url = "$repo/$asset_name"
Write-Output "Downloading: $asset_url"
New-Item -ItemType Directory -Force -Path "bqpd" | Out-Null
Set-Location "bqpd"
Invoke-WebRequest -Uri $asset_url -OutFile "$asset_name"
tar -xzvf "$asset_name"
Pop-Location

# download HiGHS
$version = "v1.11.0"
$repo = "https://github.com/amontoison/HiGHS_static_jll.jl/releases/download/HiGHS_static-${VERSION}%2B0"
$asset_name = "HiGHS_static.${VERSION}.${ARCH}-${OS}-libgfortran5.tar.gz"
$asset_url = "$repo/$asset_name"
Write-Output "Downloading: $asset_url"
New-Item -ItemType Directory -Force -Path "highs" | Out-Null
Set-Location "highs"
Invoke-WebRequest -Uri $asset_url -OutFile "$asset_name"
tar -xzvf "$asset_name"
Pop-Location