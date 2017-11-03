#!/bin/bash

for ((c=1; c <= $1; c++))
do 
echo qsub run$c.pbs >> qsub_list.sh
done
