#!/bin/bash
for file in `ls -1 fulldb/*faa`; do
	strain=`echo -e "$file" | sed 's/fulldb\///g'`
	curcommand="bash ./scripts/filterblast_ind.sh $strain"
	SGE_Batch -c "$curcommand" -P 1 -q bpp -r sge.$strain\.finaltable
done
