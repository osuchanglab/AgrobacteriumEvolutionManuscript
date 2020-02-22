#!/bin/bash
strainA=$1
declare -a seqlens
while read line;do
	c1=`echo -e "$line" | cut -f 1`
	c2=`echo -e "$line" | cut -f 2`
	seclens[c1]=$c2;
done < ./fulldb/$strain\.fai

for file in `ls -1 blast/$strainA*tab`; do
	oldid="nothing"
	hits=0
	strainB=`echo $file | sed 's/^.*_vs_//g;s/.tab//g;s/blast\///g'`
	
	while read line; do
		curid=`echo -e "$line" | cut -f 1`
		cureval=`echo -e "$line" | cut -f 11 | sed 's/e/E/g'`
		curpctid=`echo -e "$line" | cut -f 3 | sed 's/\..*$//g'`
		if [ $curid != $oldid ] && (( $(echo "$cureval < 0.00001" | bc -l) )) && [ $curpctid -gt 49 ]; then 
			querylength=${seclens[curid]} #convert to length
			curalnlen=`echo -e "$line" | cut -f 4`
			curquercov=`echo "$curalnlen/$querylength * 100" | bc -l | sed 's/\..*$//g'`
			if [ $curquercov -gt 39  ]; then
				hits=$((hits+1))
			fi
		fi
		oldid=$curid
	done < $file
	echo "$strainA	$strainB	$hits" >> pocpnew.tmp
done


