#!/usr/bin/env bash

# find location of this script
ROOT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"
TEST_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# include this project in HOUDINI_PATH
export HOUDINI_PATH="${ROOT_DIR}:&:${HOUDINI_PATH}"
#export HOUDINI_OCL_DEVICETYPE="CPU"

# source H env
cd /opt/hfs16.0.736/
source houdini_setup
cd $TEST_DIR

# run the scene
houdini testing_scene.hipnc -foreground
