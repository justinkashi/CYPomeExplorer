#!/bin/bash

# Define variables
SYS="eq_cyp"
BIN="/mms/common/gromos_20210129/gromosXX/BUILD_CUDA/bin/md"
DIR=$(pwd)
TOPO="../../topo/2HI4.top"
INPUT="../equilibration.imd"
COORD="../2HI4.cnf"
POSRESSPEC="../2HI4.pos"
REFPOS="../2HI4.rpr"
TEMPLATE="/home/teakuvek/libraries/mk_script_cuda_8_slurm-v2.lib"
VERSION="md++"
JOBLIST="../equilibration.jobs"

# Call the mk_script with provided arguments
mk_script @sys $SYS @bin $BIN @dir $DIR @files topo $TOPO input $INPUT coord $COORD posresspec $POSRESSPEC refpos $REFPOS @template $TEMPLATE @version $VERSION @joblist $JOBLIST

