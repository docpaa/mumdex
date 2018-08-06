#! /bin/bash

#
# common commands to make shell scripting better
#

set -o pipefail
function warn() {
    echo "$@" 1>&2
}
function echoe() {
    echo "$@" 1>&2
}
function echos() {
    warn "$@"
    echo "$@"
}
function report() {
    warn "$@"
    echo "$@"
}
function error() {
    warn Error: "$@"
    exit 1
}
function errors() {
    warnings Error: "$@"
    exit 1
}





