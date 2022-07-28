#!/bin/sh
# translate `complex` to `_Complex`
SCRIPT_DIR=$(dirname "$(readlink -f $0)")
[ $# -eq 1 ] || exit
ELPA_INCLUDE_DIR=$1
ELPA_GENERATED_H="${ELPA_INCLUDE_DIR}"/elpa/elpa_generated.h
MY_ELPA_GENERATED_H="${SCRIPT_DIR}"/my_elpa_generated.h

[ -f "$ELPA_GENERATED_H" ] || exit
sed 's/double complex/double _Complex/g' "$ELPA_GENERATED_H" > "$MY_ELPA_GENERATED_H"
