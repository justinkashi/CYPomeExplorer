#!/bin/bash

# Define variables
SYS="eq_cyp"
BIN="/pool/bbraun/resources/gromosXX/gromosXX/BUILD_cluster_cuda/bin/md"
DIR=$(pwd)
TOPO="../../topo/CYP90C1.top"
INPUT="equilibration.imd"
COORD="../cyp_ion.cnf"
POSRESSPEC="../cyp_ion.por"
REFPOS="../cyp_ion.rpr"
TEMPLATE="../mk_script_cuda_8_slurm.lib"
VERSION="md++"
JOBLIST="../equilibration.jobs"

# Call the mk_script with provided arguments
mk_script @sys $SYS @bin $BIN @dir $DIR @files topo $TOPO input $INPUT coord $COORD posresspec $POSRESSPEC refpos $REFPOS @template $TEMPLATE @version $VERSION @joblist $JOBLIST

