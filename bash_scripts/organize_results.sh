#!/bin/bash

################################################################################
# arrange the outputs of affinity test
#
# Dependencies:
#
# Example Usage:
#    ./organize_results.sh
#   make sure this script is put in the same folder as the affinity tests
################################################################################

set -x
cluster=peloton-ii-32tasks-core-openmpi-4.0.1
mkdir results/spherical_shell_expensive_solver/${cluster}
for i in `ls tmp/$cluster | grep output`
do
echo=1
j=`echo $i | cut -d _ -f 2-4`
cp tmp/$cluster/$i/log.txt results/spherical_shell_expensive_solver/${cluster}/output_$j  # haoyuan: here it copies the log.txt file which contains time info
done
