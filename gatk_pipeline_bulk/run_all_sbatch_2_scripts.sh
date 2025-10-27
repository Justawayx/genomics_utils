for script in $(ls scripts*/sbatch*)
do
	echo $script
	sbatch $script
done
