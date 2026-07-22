# shellcheck shell=bash
# Locate conda-forge flang_rt import libs + runtime DLLs for MSVC link and
# LoadLibrary of Fortran pot plugins. Source from Windows GHA bash steps:
#   source scripts/ci/win_flang_rt.sh
#   eon_export_flang_rt "${CONDA_PREFIX:-$PIXI_PRE}"
#
# Sets FLANG_RT_DIR and prepends it to LIB (link) and PATH (runtime DLLs).
# Exits 1 if required and the dir cannot be found.

eon_win_to_unix_path() {
  # flang -print-resource-dir often emits D:\...; Git Bash [ -d ] needs /d/...
  local p="${1//\\//}"
  if [[ "$p" =~ ^([A-Za-z]):(.*)$ ]]; then
    local drive="${BASH_REMATCH[1]}"
    local rest="${BASH_REMATCH[2]}"
    # lower-case drive for Git Bash
    drive="$(echo "$drive" | tr 'A-Z' 'a-z')"
    echo "/${drive}${rest}"
  else
    echo "$p"
  fi
}

eon_flang_rt_has_artifacts() {
  local d="$1"
  [ -d "$d" ] || return 1
  # Import libs (MSVC link) and/or runtime DLLs (LoadLibrary of pot plugins)
  local f
  for f in \
    flang_rt.runtime.dynamic.lib \
    FortranRuntime.dynamic.lib \
    flang_rt.dynamic.lib \
    flang_rt.runtime.dynamic.dll \
    FortranRuntime.dynamic.dll \
    flang_rt.dynamic.dll
  do
    if [ -f "$d/$f" ]; then
      return 0
    fi
  done
  # Some layouts only ship versioned names
  ls "$d"/flang_rt*.lib "$d"/FortranRuntime*.lib \
     "$d"/flang_rt*.dll "$d"/FortranRuntime*.dll \
     >/dev/null 2>&1
}

eon_find_flang_rt_dir() {
  # Optional roots: conda/pixi prefixes to search when -print-resource-dir fails
  local root res candidate found
  local roots=()
  for root in "$@"; do
    [ -n "$root" ] || continue
    roots+=("$(eon_win_to_unix_path "$root")")
  done

  res=""
  if command -v flang-new >/dev/null 2>&1; then
    res="$(flang-new -print-resource-dir 2>/dev/null || true)"
  fi
  if [ -z "$res" ] && command -v flang >/dev/null 2>&1; then
    res="$(flang -print-resource-dir 2>/dev/null || true)"
  fi
  if [ -n "$res" ]; then
    res="$(eon_win_to_unix_path "$res")"
    echo "flang -print-resource-dir -> $res" >&2
  fi

  local candidates=()
  if [ -n "$res" ]; then
    candidates+=(
      "${res}/lib/x86_64-pc-windows-msvc"
      "${res}/lib"
      "${res}"
    )
  fi
  for root in "${roots[@]}"; do
    # shellcheck disable=SC2086
    for candidate in \
      "${root}/Library/lib/clang"/*/lib/x86_64-pc-windows-msvc \
      "${root}/lib/clang"/*/lib/x86_64-pc-windows-msvc \
      "${root}/Library/lib" \
      "${root}/Library/bin"
    do
      candidates+=("$candidate")
    done
  done

  for candidate in "${candidates[@]}"; do
    # Unmatched globs stay literal; skip those
    case "$candidate" in
      *\**) continue ;;
    esac
    if eon_flang_rt_has_artifacts "$candidate"; then
      echo "$candidate"
      return 0
    fi
  done

  # Last resort: find import lib under roots
  for root in "${roots[@]}"; do
    [ -d "$root" ] || continue
    found="$(find "$root" \( \
      -name 'flang_rt.runtime.dynamic.lib' -o \
      -name 'FortranRuntime.dynamic.lib' -o \
      -name 'flang_rt.dynamic.lib' \
    \) 2>/dev/null | head -1 || true)"
    if [ -n "$found" ]; then
      dirname "$found"
      return 0
    fi
  done
  return 1
}

eon_export_flang_rt() {
  # Usage: eon_export_flang_rt [prefix ...] [--required]
  local required=0
  local roots=()
  local arg
  for arg in "$@"; do
    if [ "$arg" = "--required" ]; then
      required=1
    elif [ -n "$arg" ]; then
      roots+=("$arg")
    fi
  done

  local dir
  if ! dir="$(eon_find_flang_rt_dir "${roots[@]}")"; then
    echo "WARNING: flang_rt resource dir not found (roots: ${roots[*]:-none})" >&2
    if [ -n "${CONDA_PREFIX:-}" ]; then
      echo "CONDA_PREFIX=$CONDA_PREFIX" >&2
      ls -la "${CONDA_PREFIX}/Library/lib/clang" 2>/dev/null || true
    fi
    if command -v flang-new >/dev/null 2>&1; then
      echo "flang-new=$(command -v flang-new)" >&2
      flang-new -print-resource-dir 2>&1 || true
    fi
    if [ "$required" -eq 1 ]; then
      return 1
    fi
    return 0
  fi

  export FLANG_RT_DIR="$dir"
  # LIB uses Windows-style ';' separators for MSVC link.exe
  export LIB="${FLANG_RT_DIR};${LIB:-}"
  # Runtime DLLs for LoadLibrary of eon_*.dll pot plugins
  export PATH="${FLANG_RT_DIR}:${PATH}"
  echo "Using flang_rt LIBPATH+PATH: $FLANG_RT_DIR"
  ls "$FLANG_RT_DIR"/*.{lib,dll} 2>/dev/null | head -30 || ls "$FLANG_RT_DIR" | head -30
  return 0
}
