# Activate MSVC for GitHub Actions without clobbering the pixi/conda env.
# Clears conda-forge vs* leftovers that break VsDevCmd when vswhere selects VS 18.

$kill = @(
  'VSINSTALLDIR','VS_VERSION','VS_MAJOR','VS_YEAR','VCVARSBAT',
  'NEWER_VS_WITH_OLDER_VC','WindowsSdkDir','WindowsSDKVer','MSSdk',
  'DISTUTILS_USE_SDK','VCINSTALLDIR','VCToolsInstallDir','VCToolsVersion',
  'CMAKE_GENERATOR','CMAKE_GENERATOR_TOOLSET','CMAKE_GENERATOR_PLATFORM',
  'CMAKE_GEN','CMAKE_PLAT','USE_NEW_CMAKE_GEN_SYNTAX','VSCMD_VER',
  'VSCMD_ARG_TGT_ARCH','VSCMD_ARG_HOST_ARCH','UniversalCRTSdkDir',
  'UCRTVersion','WindowsLibPath','WindowsSdkLibPath','WindowsSdkVerBinPath',
  'WindowsSDKVersion','ExtensionSdkDir','VCToolsRedistDir','FrameworkDir',
  'FrameworkDir64','FrameworkVersion','FrameworkVersion64','Framework40Version'
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

# Snapshot env before/after vcvars; only push MSVC-related deltas + PATH/INCLUDE/LIB*
$before = @{}
cmd.exe /c set | ForEach-Object {
  if ($_ -match '^([^=]+)=(.*)$') { $before[$Matches[1]] = $Matches[2] }
}

$afterLines = cmd.exe /c "`"$vcvars`" >nul 2>&1 && set"
$after = @{}
foreach ($line in $afterLines) {
  if ($line -match '^([^=]+)=(.*)$') { $after[$Matches[1]] = $Matches[2] }
}

# Always take these from vcvars (required for cl/link)
$must = @('PATH','INCLUDE','LIB','LIBPATH','VCINSTALLDIR','VCToolsInstallDir',
  'VCToolsVersion','WindowsSdkDir','WindowsSDKVersion','WindowsSdkVerBinPath',
  'UniversalCRTSdkDir','UCRTVersion','WindowsLibPath','WindowsSdkLibPath',
  'VSINSTALLDIR','VSCMD_VER','VSCMD_ARG_TGT_ARCH','VSCMD_ARG_HOST_ARCH',
  'FrameworkDir','FrameworkDir64','FrameworkVersion','FrameworkVersion64',
  'DevEnvDir','ExtensionSdkDir','VCToolsRedistDir','WindowsSdkBinPath',
  'WindowsSdkDir','IFCPATH','NETFXSDKDir','ExternalCompilerOptions')

foreach ($name in $must) {
  if ($after.ContainsKey($name)) {
    $value = $after[$name]
    Add-Content -Path $env:GITHUB_ENV -Value "$name=$value"
    Set-Item -Path "Env:$name" -Value $value
  }
}

# Also export any VSCMD_* / VC* / WindowsSDK* newly set by vcvars
foreach ($name in $after.Keys) {
  if ($name -match '^(VSCMD_|VC|WindowsSDK|WindowsSdk|UniversalCRT|Framework|DevEnv|IFC|NETFX)') {
    if (-not ($must -contains $name)) {
      $value = $after[$name]
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
