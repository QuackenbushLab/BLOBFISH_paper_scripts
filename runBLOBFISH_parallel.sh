#!/bin/bash
# For each file in the expression directory, run Panda and save results.
# Create a function that runs the next process only when memory is available.
randoutputs=""
check_memory_and_run() {
	available_memory=$(free -m | awk '/^Mem:/{print $4}')
	if [ ! -e "$randoutputs/$1_$2.RDS" ]; then
		#if ((available_memory > 20000)); then
			echo "running blobfish"
			echo $1
			echo $2
			echo $3
			sudo Rscript single_blobfish_run.R $1 $2 $3 > $1$2.txt &
			return 0
		#else
		#	echo "Insufficient memory available. Waiting..."
		#	return 1
		#fi
		echo "running"
	else echo "$randoutputs/$1_$2.RDS"
	fi
}

# For each file, check if memory is available, then run.
tissues=("adipose" "aorta" "lung" "muscle" "skin")
pvalsdir=""
pvals=("$pvalsdir/adiposePval.RDS" "$pvalsdir/aortaPval.RDS" "$pvalsdir/lungPval.RDS" "$pvalsdir/skeletalMusclePval.RDS" "$pvalsdir/skinPval.RDS")
pvals=("/home/ubuntu/BLOBFISH/skinPval.RDS")
for (( i=0;i<${#tissues[@]}; i++ ));
do
	for (( j=1;j<100; j++ ));
	do
		echo ${tissues[$i]}
		while ! check_memory_and_run ${tissues[$i]} $j ${pvals[$i]};
		do
			sleep 60
		done
		sleep 60
	done
done

# Wait until all jobs have completed.
wait
echo "All networks ready!"

