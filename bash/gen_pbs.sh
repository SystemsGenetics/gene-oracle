#!/bin/bash

for ((c=1; c <= $1; c++))
do
var2="$c * 256"
var=$(($var2))
cat > run$c.pbs << EOF$c 
#PBS -N hidden_layer_test
#PBS -l select=1:ncpus=16:ngpus=2:mem=16gb:gpu_model=k40,walltime=24:00:00
#PBS -M efwoods@clemson.edu
#PBS -m a
cd /home/efwoods/DeepGTEx                                                       

singularity exec --nv -B /software:/software -B /scratch2:/scratch2 -B /local_scratch:/local_scratch /software/singularity-containers/ubuntu_dl_gpu.img /usr/bin/python ~/DeepGTEx/nn_gtex.py --epochs 150 --h1 $var --h2 $var --h3 $var --h4 $var --n_classes 53  > ~/DeepGTEx/logs/gtex_$varx$varx$varx$var_lr_001_e_150_$c
EOF$c
done
