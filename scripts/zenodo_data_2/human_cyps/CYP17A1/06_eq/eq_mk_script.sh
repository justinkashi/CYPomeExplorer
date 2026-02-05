#!/bin/bash

# Define variables
SYS="eq_cyp"
BIN="/mms/common/gromos_20231119/md++/BUILD_deb12_cuda/bin/md"
DIR=$(pwd)
TOPO="../../topo/CYP17A1.top"
INPUT="equilibration.imd"
COORD="../cyp_ion.cnf"
POSRESSPEC="../cyp_ion.por"
REFPOS="../cyp_ion.rpr"
TEMPLATE="/home/teakuvek/libraries/mk_script_cuda_8_slurm-v2.lib"
VERSION="md++"
JOBLIST="../equilibration.jobs"

# Call the mk_script with provided arguments
mk_script @sys $SYS @bin $BIN @dir $DIR @files topo $TOPO input $INPUT coord $COORD posresspec $POSRESSPEC refpos $REFPOS @template $TEMPLATE @version $VERSION @joblist $JOBLIST

