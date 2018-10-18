#!/usr/bin/env bash

# usage

# find location of this script and move to blink dir
ROOT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd ../blink && pwd )"

# set needed paths
export BLINK_PATH="${ROOT_DIR}"
export FN_BLINK_INCLUDE_PATHS="${ROOT_DIR}"

# run the scene, expects an "nuke_launcher" environment variable pointing to binary of Nuke (e.g. /usr/local/Nuke11.1v4/Nuke11.1)
${nuke_launcher} test_scene.nk
