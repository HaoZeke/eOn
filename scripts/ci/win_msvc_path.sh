# shellcheck shell=bash
# Source from pixi-run bash steps on Windows GHA.
# Prefer conda DLL dirs (torch), then MSVC cl/link, never Git /usr/bin link.
_path_strip_usr_bin() {
  local _out="" _p
  local IFS=':'
  for _p in $PATH; do
    [ "$_p" = "/usr/bin" ] && continue
    if [ -z "$_out" ]; then
      _out="$_p"
    else
      _out="$_out:$_p"
    fi
  done
  PATH="$_out"
  export PATH
}

_path_strip_usr_bin
export PATH="${CONDA_PREFIX}/Scripts:${CONDA_PREFIX}:${CONDA_PREFIX}/Library/bin:${PATH}"
_cl="$(command -v cl.exe 2>/dev/null || true)"
if [ -n "$_cl" ]; then
  _cldir="$(dirname "$_cl")"
  if [ -x "${_cldir}/link.exe" ]; then
    export PATH="${_cldir}:${PATH}"
  fi
fi
unset LDFLAGS || true
export LD=link
export CC=cl
export CXX=cl
