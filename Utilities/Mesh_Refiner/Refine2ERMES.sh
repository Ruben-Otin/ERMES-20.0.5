#!/bin/bash
#==============================================================================
# Refine2ERMES.sh - refine an ERMES mesh and run ERMES on the refined mesh
#==============================================================================
#
# Pipeline (all inside the GiD problem folder <name>.gid):
#   1. remove the results of the previous run
#   2. convert the .dat files written by GiD to Unix line endings (dos2unix)
#   3. refine the mesh with mesh_refiner*.py  (rewrites <name>-N.dat in place)
#      (skipped when the mesh is already refined, or with FORCE_REFINE=2)
#   4. run ERMES on the refined mesh          (log -> <name>.info)
#
# Usage (from a WSL terminal, or from the GiD problem-type .bat through wsl):
#   ./Refine2ERMES.sh                       # defaults of the USER SETTINGS below
#   ./Refine2ERMES.sh MyProblem             # other problem name
#   ./Refine2ERMES.sh MyProblem 3           # ... refined with LEVEL = 3
#   LEVEL=4 NPROC=8 ./Refine2ERMES.sh       # any setting can be given as an
#                                           # environment variable instead
#   REFINER=/path/mesh_refiner.py ERMES_EXE=/path/ERMES_x ./Refine2ERMES.sh
#                                           # problem, refiner and ERMES can be
#                                           # in three different locations
#   FORCE_REFINE=2 ./Refine2ERMES.sh        # run ERMES alone, never refine
#   ./Refine2ERMES.sh -h                    # this help
#
# IMPORTANT - the refinement is NOT idempotent. mesh_refiner rewrites the .dat
# files in place, so running it twice would refine the already refined mesh
# (tetrahedra x LEVEL^3 each time). To prevent that, a stamp file
# (.refined_stamp) is created after a successful refinement; while no .dat file
# is newer than the stamp the refinement is skipped and ERMES simply runs again.
# Every GiD "Calculate" rewrites the .dat files, which re-enables the refinement.
#
# FORCE_REFINE selects what happens with the refinement:
#   0 : default - automatic: refine only if a .dat file is newer than the stamp
#   1 : always refine, even if the mesh looks already refined
#   2 : never refine: ERMES is run alone on the .dat files exactly as
#       they are now, even if they have changed (e.g. a refined mesh
#       edited by hand). The stamp is left untouched, so a later run
#       with FORCE_REFINE=0 will refine again if the .dat files are newer.
#==============================================================================

# Stop at the first error, on unset variables and on failures inside pipes.
# Without this a failed refinement would still launch ERMES on the old mesh.
set -euo pipefail
trap 'echo "ERROR: $(basename "$0") failed at line $LINENO (exit code $?)" >&2' ERR

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    # print the header comment (everything up to the first non-comment line)
    awk 'NR>1 && /^#/ {sub(/^# ?/, ""); print; next} NR>1 {exit}' "$0"
    exit 0
fi

#==============================================================================
# USER SETTINGS
#==============================================================================
# Every setting is resolved in this order of priority:
#     command-line argument  >  environment variable  >  DEF_* default below
# Edit only the DEF_* defaults of the DEFAULT VALUES block.

# ---------------------------- DEFAULT VALUES ----------------------------------

# Paths are WSL paths (D:\ON-WORK\... = /mnt/d/ON-WORK/...). Linux paths are case
# sensitive. The three locations are independent of each other.

DEF_PROBLEM_NAME=MyProblemName               # GiD project name without ".gid"
DEF_PROBLEM_FOLDER=/path_to_work/Folder      # folder that CONTAINS <name>.gid

DEF_REFINER=/path_to/mesh_refiner.py         # full path of the refiner script
DEF_ERMES_EXE=/path_to/ERMES_20.0.5.exe      # full path of the ERMES executable

DEF_MODE=levels                              # levels | ps
DEF_LEVEL=2                                  # edge divisions (<= 1: no refinement)
DEF_NPROC=$(nproc)                           # refiner processes (default: all cores)
DEF_TMP_DIR=/tmp                             # temporary files of the refiner
DEF_FORCE_REFINE=0                           # 0=auto | 1=force refine | 2=only ERMES

# ------------------------- RESOLVED SETTINGS ----------------------------------

# Problem name = GiD project name without ".gid"  (argument 1)
PROBLEM_NAME=${1:-${PROBLEM_NAME:-$DEF_PROBLEM_NAME}}

# Folder that contains the GiD project folder(s) <name>.gid, and the project folder
# actually used. The latter follows the PROBLEM_NAME in use (argument 1 / environment).
PROBLEM_FOLDER=${PROBLEM_FOLDER:-$DEF_PROBLEM_FOLDER}
PROBLEM_FOLDER=${PROBLEM_FOLDER%/}               # drop a trailing slash, if any
PROBLEM_DIR=$PROBLEM_FOLDER/$PROBLEM_NAME.gid

# Refiner script and ERMES executable: full paths, in any location
REFINER=${REFINER:-$DEF_REFINER}
ERMES_EXE=${ERMES_EXE:-$DEF_ERMES_EXE}

# Refinement: MODE = levels (edges divided in LEVEL parts: tets x LEVEL^3) or
#                    ps     (Powell-Sabin: tets x 24, LEVEL is ignored)
# LEVEL <= 1 -> no refinement, ERMES is run on the mesh exactly as GiD wrote it.
MODE=${MODE:-$DEF_MODE}
LEVEL=${2:-${LEVEL:-$DEF_LEVEL}}

# Parallel processes for the refiner
NPROC=${NPROC:-$DEF_NPROC}

# Folder for the temporary part files of the refiner. A Linux folder is MUCH
# faster than the Windows drive (/mnt/d). It needs ~ the size of the refined .dat.
TMP_DIR=${TMP_DIR:-$DEF_TMP_DIR}

# 0 = refine only if needed | 1 = always refine | 2 = never refine, ERMES alone
# (see IMPORTANT above)
FORCE_REFINE=${FORCE_REFINE:-$DEF_FORCE_REFINE}

#==============================================================================
# CHECKS
#==============================================================================

# numeric settings must be plain integers (otherwise bash arithmetic fails cryptically)
for v in LEVEL NPROC; do
    [[ "${!v}" =~ ^[0-9]+$ ]] || { echo "ERROR: $v must be a non-negative integer (got '${!v}')" >&2; exit 1; }
done
[[ "$FORCE_REFINE" =~ ^[012]$ ]] || { echo "ERROR: FORCE_REFINE must be 0, 1 or 2 (got '$FORCE_REFINE')" >&2; exit 1; }
[[ "$MODE" == "levels" || "$MODE" == "ps" ]] || { echo "ERROR: MODE must be 'levels' or 'ps' (got '$MODE')" >&2; exit 1; }
(( NPROC >= 1 )) || { echo "ERROR: NPROC must be >= 1" >&2; exit 1; }

[[ -d "$PROBLEM_DIR" ]] || { echo "ERROR: problem folder not found: $PROBLEM_DIR" >&2; exit 1; }
PROBLEM_DIR=$(cd "$PROBLEM_DIR" && pwd)          # absolute path (the script cd's below)

# Refiner and ERMES may be relative paths: make them absolute now, because the
# script changes to the problem folder before using them.
REFINER_PATH=$(realpath -m "$REFINER")

if [[ "$ERMES_EXE" == */* ]]; then
    ERMES_PATH=$(realpath -m "$ERMES_EXE")
else                                             # bare command name: look it up in PATH
    ERMES_PATH=$(command -v "$ERMES_EXE" || true)
fi
[[ -n "$ERMES_PATH" && -x "$ERMES_PATH" && ! -d "$ERMES_PATH" ]] || \
    { echo "ERROR: ERMES executable not found or not executable: $ERMES_EXE" >&2; exit 1; }

# the refiner (and python3) are only needed when a refinement can take place
if [[ "$FORCE_REFINE" != "2" ]] && { (( LEVEL > 1 )) || [[ "$MODE" == "ps" ]]; }; then
    [[ -f "$REFINER_PATH" ]] || { echo "ERROR: refiner not found: $REFINER_PATH" >&2; exit 1; }
    command -v python3 >/dev/null || { echo "ERROR: python3 not found" >&2; exit 1; }
fi

cd "$PROBLEM_DIR"
echo "Problem : $PROBLEM_NAME  ($PROBLEM_DIR)"

# nullglob: a pattern without matches expands to nothing (instead of itself)
shopt -s nullglob
DAT_FILES=("$PROBLEM_NAME"*.dat)
(( ${#DAT_FILES[@]} > 0 )) || { echo "ERROR: no $PROBLEM_NAME*.dat files in $PROBLEM_DIR (run Calculate in GiD first)" >&2; exit 1; }

#==============================================================================
# 1. REMOVE OLD RESULTS
#==============================================================================

# (-f: no error if a file does not exist)
rm -f "$PROBLEM_NAME.info" \
      "$PROBLEM_NAME.flavia.res" \
      "$PROBLEM_NAME.post.res" \
      "$PROBLEM_NAME.post.msh"

#==============================================================================
# 2. LINE ENDINGS  (always, also when the refinement is skipped)
#==============================================================================

# GiD (Windows) writes CRLF line endings; convert to LF for the Linux tools.
# Only the files that really contain a CR are touched, and their modification
# time is preserved (dos2unix -k, or touch -r with the sed fallback) so that the
# ".refined_stamp" comparison below is not disturbed.
for f in "${DAT_FILES[@]}"; do
    if grep -q $'\r' "$f"; then
        if command -v dos2unix >/dev/null; then
            dos2unix -q -k "$f"
        else
            ref=$(mktemp); touch -r "$f" "$ref"
            sed -i 's/\r$//' "$f"
            touch -r "$ref" "$f"; rm -f "$ref"
        fi
    fi
done

#==============================================================================
# 3. MESH REFINEMENT  (skipped if already refined, or with FORCE_REFINE=2)
#==============================================================================

STAMP=.refined_stamp
DO_REFINE=1
if [[ "$FORCE_REFINE" == "2" ]]; then
    echo "Refinement disabled (FORCE_REFINE=2): ERMES runs on the .dat files as they are."
    DO_REFINE=0
elif (( LEVEL <= 1 )) && [[ "$MODE" != "ps" ]]; then
    echo "Refinement skipped (LEVEL = $LEVEL)"
    DO_REFINE=0
elif [[ -f "$STAMP" && "$FORCE_REFINE" == "0" ]] \
     && [[ -z "$(find . -maxdepth 1 -name "$PROBLEM_NAME*.dat" -newer "$STAMP" -print -quit)" ]]; then
    echo "Mesh already refined (no .dat file newer than $STAMP): refinement skipped."
    echo "  -> run Calculate in GiD again, or use FORCE_REFINE=1, to refine again."
    DO_REFINE=0
fi

if (( DO_REFINE )); then
    # Refine. --name avoids the auto-detection of the problem name,
    #    --no-backup because GiD regenerates the originals at every Calculate.
    REFINE_ARGS=("$PROBLEM_DIR" --name "$PROBLEM_NAME" --mode "$MODE"
                 --np "$NPROC" --tmp "$TMP_DIR" --no-backup)
    [[ "$MODE" == "levels" ]] && REFINE_ARGS+=(--level "$LEVEL")

    echo "Refining: python3 $(basename "$REFINER_PATH") ${REFINE_ARGS[*]}"
    # tee keeps a copy of the screen output (pipefail makes a refiner failure abort)
    python3 "$REFINER_PATH" "${REFINE_ARGS[@]}" 2>&1 | tee "$PROBLEM_NAME.refine.log"

    touch "$STAMP"          # refined: remember it (only reached if the refiner succeeded)
fi

#==============================================================================
# 4. RUN ERMES
#==============================================================================

# stdout and stderr both go to <name>.info
echo "Running ERMES: $ERMES_PATH $PROBLEM_NAME   (log: $PROBLEM_NAME.info)"
SECONDS=0
"$ERMES_PATH" "$PROBLEM_NAME" > "$PROBLEM_DIR/$PROBLEM_NAME.info" 2>&1
echo "ERMES finished in ${SECONDS} s"
