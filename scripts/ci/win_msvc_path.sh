# shellcheck shell=bash
# Windows GHA: make MSVC available for meson without breaking torch DLL load.
#
# Order matters:
# 1) Capture absolute cl/link while vcvars PATH is still present
# 2) Rebuild PATH with conda + torch/lib first (shm.dll deps)
# 3) Append MSVC tool dir after conda so Python loads conda CRTs first
# 4) Drop Git /usr/bin so GNU link.exe cannot win

_path_without_usr_bin() {
  local _out="" _p
  local IFS=':'
  for _p in $1; do
    [ "$_p" = "/usr/bin" ] && continue
    [ -z "$_p" ] && continue
    if [ -z "$_out" ]; then
      _out="$_p"
    else
      _out="${_out}:${_p}"
    fi
  done
  printf '%s' "$_out"
}

# Capture MSVC tools before we rewrite PATH
_CL="$(command -v cl.exe 2>/dev/null || true)"
_LINK="$(command -v link.exe 2>/dev/null || true)"
# Prefer real MSVC link if GNU link is first
if [ -n "$_CL" ]; then
  _CLDIR="$(dirname "$_CL")"
  if [ -x "${_CLDIR}/link.exe" ]; then
    _LINK="${_CLDIR}/link.exe"
  fi
fi

PATH="$(_path_without_usr_bin "$PATH")"
export PATH

# Conda / pixi first (DLLs for torch, numpy, etc.)
if [ -n "${CONDA_PREFIX:-}" ]; then
  export PATH="${CONDA_PREFIX}/Scripts:${CONDA_PREFIX}:${CONDA_PREFIX}/Library/bin:${PATH}"
  # torch native libs (shm.dll and friends)
  _torch_lib="${CONDA_PREFIX}/Lib/site-packages/torch/lib"
  if [ -d "$_torch_lib" ]; then
    export PATH="${_torch_lib}:${PATH}"
  fi
fi

# MSVC tool dir after conda so LoadLibrary prefers conda CRT/OpenMP
if [ -n "$_CL" ]; then
  export CC="$_CL"
  export CXX="$_CL"
  export PATH="$(dirname "$_CL"):${PATH}"
fi
if [ -n "$_LINK" ]; then
  export LD="$_LINK"
else
  export LD=link
fi

# Drop conda clang/lld flags that break MSVC link lines
unset LDFLAGS || true
# Do not force lld via FFLAGS for C++ link
if [ -n "${FFLAGS:-}" ]; then
  FFLAGS="$(printf '%s' "$FFLAGS" | sed 's/-fuse-ld=lld//g')"
  export FFLAGS
fi

# Sanity
command -v cl.exe >/dev/null
command -v link.exe >/dev/null
# link must not be GNU
if link --help 2>&1 | head -1 | grep -qi 'Usage: link'; then
  : # MSVC link --help says "usage: LINK" sometimes; check differently
  true
fi
if command -v link.exe | grep -Eqi 'Git|usr/bin|mingw'; then
  echo "ERROR: link.exe is still GNU: $(command -v link.exe)" >&2
  exit 1
fi
