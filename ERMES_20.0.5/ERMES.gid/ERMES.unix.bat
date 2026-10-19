#!/bin/bash

# ProblemName     = $1
# ProblemPath     = $2
# ProblemTypePath = $3

# OutputFile: $2/$1.info

# ERMES executable. To use another build, replace this file or change the name here.
ERMES_EXE="$3/execs/unix/ERMES_20.0.5"

rm -f "$2/$1.info" 
rm -f "$2/$1.flavia.res" 
rm -f "$2/$1.post.res"
rm -f "$2/$1.post.msh"

# Problemtypes unpacked from a zip file lose the execute permission of the executable
[ -x "$ERMES_EXE" ] || chmod +x "$ERMES_EXE" 2>/dev/null

# Run from the problem folder ( ERMES reads and writes its files in the working directory ).
# Error messages ( stderr ) are also written to the .info file shown by GiD.
cd "$2" || exit 1

"$ERMES_EXE" "$1" > "$2/$1.info" 2>&1
