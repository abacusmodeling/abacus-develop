#!/bin/sh
SCRIPT_DIR=$(dirname "$(readlink -f $0)")
[ $# -eq 1 ] || exit
ELPA_INCLUDE_DIR=$1
ELPA_GENERIC_H="${ELPA_INCLUDE_DIR}"/elpa/elpa_generic.h
MY_ELPA_GENERIC_HPP="${SCRIPT_DIR}"/my_elpa_generic.hpp

[ -f "$ELPA_GENERIC_H" ] || exit
# check whether the file `elpa_generic.h` has keywords `elpa_eigenvectors_all_host_arrays_dc`
# if it has, it is the new version in 2021.11.002; otherwise it is the old version
if [ "$(grep elpa_eigenvectors_all_host_arrays_dc $ELPA_GENERIC_H)" ]
then
    cp "${SCRIPT_DIR}"/elpa_generic_template_2.hpp "$MY_ELPA_GENERIC_HPP"
else
    cp "${SCRIPT_DIR}"/elpa_generic_template_1.hpp "$MY_ELPA_GENERIC_HPP"
fi
