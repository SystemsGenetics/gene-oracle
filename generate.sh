~/DeepGTEx/gen_pbs.sh $1
echo > ~/DeepGTEx/qsub_list.sh
~/DeepGTEx/sub_pbs.sh $1
~/DeepGTEx/qsub_list.sh
