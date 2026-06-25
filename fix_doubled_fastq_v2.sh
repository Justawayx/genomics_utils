#!/bin/bash

for filename in $(ls *fastq.gz)
do
	echo "Working on ${filename}..."
	NUM_LINES=$(zcat $filename | wc -l)
	echo "${filename} has ${NUM_LINES} lines"
	zcat ${filename} | head -n $((${NUM_LINES} / 2)) | pigz > ${filename}_fixed
	echo "Done!"
done
