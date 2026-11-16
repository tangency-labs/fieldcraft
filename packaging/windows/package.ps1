<#
.SYNOPSIS
  Build and package FieldCraft for Windows x64.

.DESCRIPTION
  Produces, in $env:DIST (default: dist/release):
    fieldcraft-<version>-windows-x64.zip   fieldcraft.exe, fieldcraft-cli.exe, README.md and the
                                           licence files (LICENSE-MIT, LICENSE-APACHE, NOTICE,
                                           THIRD-PARTY-LICENSES.txt)

  The binaries link the C runtime statically (+crt-static, .cargo/config.toml), so nothing else
  needs installing.
  fieldcraft.exe carries the application icon (assets/app-icon/icon.ico, embedded by
  apps/fieldcraft/build.rs). Nothing is signed.

  Adapted from SolidCraft's packaging/windows/package.ps1 (MIT OR Apache-2.0; see ATTRIBUTION.md).

.EXAMPLE
  pwsh packaging/windows/package.ps1
#>
param([switch] $SkipBuild)
$ErrorActionPreference = 'Stop'
$Root = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path

function Invoke-Native([string] $What, [scriptblock] $Block) {
  Write-Output "==> $What"
  & $Block
  if ($LASTEXITCODE -ne 0) { throw "$What failed with exit code $LASTEXITCODE" }
}

# The version: FIELDCRAFT_VERSION (CI sets it), else [workspace.package] version in Cargo.toml.
$Version = $env:FIELDCRAFT_VERSION
if (-not $Version) {
  $meta = cargo metadata --no-deps --format-version 1 --manifest-path (Join-Path $Root 'Cargo.toml') | ConvertFrom-Json
  $Version = ($meta.packages | Where-Object name -eq 'fieldcraft').version
}
if (-not $Version) { throw 'could not read [workspace.package] version from Cargo.toml' }

$Arch = 'x64'
$Target = 'x86_64-pc-windows-msvc'
$Dist = if ($env:DIST) { $env:DIST } else { Join-Path $Root 'dist\release' }
$TargetDir = if ($env:CARGO_TARGET_DIR) { $env:CARGO_TARGET_DIR } else { Join-Path $Root 'target' }
New-Item -ItemType Directory -Force -Path $Dist | Out-Null
Write-Output "FieldCraft $Version for Windows $Arch ($Target)"

if (-not $SkipBuild) {
  Invoke-Native "cargo build ($Target)" { cargo build --release --locked -p fieldcraft -p fieldcraft-cli --target $Target }
}
$Bin = Join-Path $TargetDir "$Target\release"

$Check = Join-Path $PSScriptRoot 'check-exe.ps1'
& $Check -Path (Join-Path $Bin 'fieldcraft.exe') -Subsystem 2 -Icon
& $Check -Path (Join-Path $Bin 'fieldcraft-cli.exe') -Subsystem 3

# ---- zip ---------------------------------------------------------------------------------------
$Name = "fieldcraft-$Version-windows-$Arch"
$Stage = Join-Path $TargetDir "windows-package\$Name"
Remove-Item -Recurse -Force $Stage -ErrorAction SilentlyContinue
New-Item -ItemType Directory -Force -Path $Stage | Out-Null
Copy-Item (Join-Path $Bin 'fieldcraft.exe'), (Join-Path $Bin 'fieldcraft-cli.exe') $Stage
# The files copy_docs (packaging/env.sh) puts in every package.
foreach ($f in 'README.md', 'LICENSE-MIT', 'LICENSE-APACHE', 'NOTICE', 'THIRD-PARTY-LICENSES.txt') { Copy-Item (Join-Path $Root $f) $Stage }
$Zip = Join-Path $Dist "$Name.zip"
Remove-Item -Force $Zip -ErrorAction SilentlyContinue
# Windows' own bsdtar writes the zip (-a: the format from the name; Compress-Archive is slow). By
# full path, as Git for Windows puts a GNU tar without zip support on the PATH.
$Tar = Join-Path $env:SystemRoot 'System32\tar.exe'
Invoke-Native 'zip' { & $Tar -a -cf $Zip -C (Split-Path $Stage) $Name }

# The CLI reports the version it was packaged as (a nightly's suffix included). The app is a GUI
# program (no console), so its --version is checked on Linux and macOS.
$Reported = & (Join-Path $Stage 'fieldcraft-cli.exe') --version
if ($LASTEXITCODE -ne 0) { throw "fieldcraft-cli --version failed with exit code $LASTEXITCODE" }
if ("$Reported".Trim() -ne "fieldcraft-cli $Version") { throw "fieldcraft-cli --version printed '$Reported', not 'fieldcraft-cli $Version'" }
Get-Item $Zip | Format-Table Name, Length
