# Clear conda-forge vs* leftovers that break VsDevCmd when the runner has VS 18+
# and vswhere selects a different install than conda's VSINSTALLDIR.
$kill = @(
  'VSINSTALLDIR','VS_VERSION','VS_MAJOR','VS_YEAR','VCVARSBAT',
  'NEWER_VS_WITH_OLDER_VC','WindowsSdkDir','WindowsSDKVer','MSSdk',
  'DISTUTILS_USE_SDK','VCINSTALLDIR','VCToolsInstallDir','VCToolsVersion',
  'CMAKE_GENERATOR','CMAKE_GENERATOR_TOOLSET','CMAKE_GENERATOR_PLATFORM',
  'CMAKE_GEN','CMAKE_PLAT','USE_NEW_CMAKE_GEN_SYNTAX','VSCMD_VER',
  'VSCMD_ARG_TGT_ARCH','VSCMD_ARG_HOST_ARCH','UniversalCRTSdkDir',
  'UCRTVersion','WindowsLibPath','WindowsSdkLibPath','WindowsSdkVerBinPath',
  'WindowsSDKVersion','ExtensionSdkDir','VCToolsRedistDir'
)
foreach ($v in $kill) {
  Remove-Item "Env:$v" -ErrorAction SilentlyContinue
}

$vswhere = Join-Path ${env:ProgramFiles(x86)} "Microsoft Visual Studio\Installer\vswhere.exe"
if (-not (Test-Path $vswhere)) { throw "vswhere not found at $vswhere" }

$installPath = & $vswhere -latest -products * `
  -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 `
  -property installationPath
if (-not $installPath) { throw "vswhere: no VC tools installation" }

$vcvars = Join-Path $installPath "VC\Auxiliary\Build\vcvars64.bat"
if (-not (Test-Path $vcvars)) { throw "vcvars64.bat missing: $vcvars" }

Write-Host "Using $vcvars"
# Capture environment after vcvars64
$cmdOut = cmd.exe /c "`"$vcvars`" >nul 2>&1 && set"
foreach ($line in $cmdOut) {
  if ($line -match '^([^=]+)=(.*)$') {
    $name = $Matches[1]
    $value = $Matches[2]
    # Skip empty names; append to GITHUB_ENV for subsequent steps
    if ($name) {
      Add-Content -Path $env:GITHUB_ENV -Value "$name=$value"
      Set-Item -Path "Env:$name" -Value $value
    }
  }
}

$cl = Get-Command cl.exe -ErrorAction SilentlyContinue
$link = Get-Command link.exe -ErrorAction SilentlyContinue
if (-not $cl) { throw "cl.exe not on PATH after vcvars64" }
if (-not $link) { throw "link.exe not on PATH after vcvars64" }
Write-Host "cl: $($cl.Source)"
Write-Host "link: $($link.Source)"
if ($link.Source -match 'Git|usr\\bin|mingw') {
  throw "link.exe is still GNU/Git, not MSVC: $($link.Source)"
}
