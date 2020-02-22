#!/bin/bash

for file in `ls -1 ../fulldb/*faa`; do
	echo $file
	samtools faidx $file
done
