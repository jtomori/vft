#!/usr/bin/env bash

# usage
printf "Usage:./run.sh scene_file.hip (-cpu)\n  -cpu: optional argument to run OpenCL code on a CPU device (default: GPU)\n\n\n"

# find location of this script
ROOT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd ../.. && pwd )"

# include this project in HOUDINI_PATH
export HOUDINI_PATH="${ROOT_DIR}/houdini:&:${HOUDINI_PATH}"

# setup OSL
export ARNOLD_PLUGIN_PATH="${ARNOLD_PLUGIN_PATH}:${ROOT_DIR}/osl:${ROOT_DIR}/houdini/ocl/include"
export OSL_OPTIONS="range_checking=1,debug_uninit=1,max_warnings_per_thread=100,compile_report=1"

# use a cpu mode
if [[ $2 == "-cpu" ]]
    then
        export HOUDINI_OCL_DEVICETYPE="CPU"
        printf "CPU mode set\n"
fi

# run the scene
houdini $1 -foreground