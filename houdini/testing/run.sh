#!/usr/bin/env bash

# usage
printf "Usage:\n  -cpu argument to run OpenCL code on a CPU device (by default will run on GPU)\n\n\n"

# find location of this script
ROOT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"
TEST_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# include this project in HOUDINI_PATH
export HOUDINI_PATH="${ROOT_DIR}:&:${HOUDINI_PATH}"

# use a cpu mode
if [[ $1 == "-cpu" ]]
    then
        export HOUDINI_OCL_DEVICETYPE="CPU"
        printf "CPU mode set\n"
fi

# source H env
cd /opt/hfs16.0.736/
source houdini_setup
cd $TEST_DIR

# run the scene
houdini testing_scene.hipnc -foreground
