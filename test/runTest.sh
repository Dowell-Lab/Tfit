#!/bin/bash

date
mpirun -np 1 ../src/Tfit model -v -config ../config_files/config_file.txt -ij ../../../data/PRO-BSA-Eric-1.bedGraph_typical_region.bg -N testmodel -k singleregion.bed -o .
date
