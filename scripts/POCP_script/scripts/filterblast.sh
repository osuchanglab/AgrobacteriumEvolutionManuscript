#!/bin/bash
for file in `ls -1 blast/*tab`; do
	oldid="nothing"
	hits=0
	#A4.faa_vs_Agrobacterium_tumefaciens_5A.faa.tab
	strainA=`echo $file | sed 's/_vs_.*//g;s/blast\///g'`
	strainB=`echo $file | sed 's/^.*_vs_//g;s/.tab//g;s/blast\///g'`
	
	declare -a seqlens
	while read line;do 
		c1=`echo -e "$line" | cut -f 1`
		c2=`echo -e "$line" | cut -f 2`
		seclens[c1]=$c2;
	done < ./fulldb/$strainA\.fai
	while read line; do
		curid=`echo -e "$line" | cut -f 1`
		cureval=`echo -e "$line" | cut -f 11 | sed 's/e/E/g'`
		curpctid=`echo -e "$line" | cut -f 3 | sed 's/\..*$//g'`
		if [ $curid != $oldid ] && (( $(echo "$cureval < 0.00001" | bc -l) )) && [ $curpctid -gt 49 ]; then 
			#echo -e "$line" >> $file\.filtered.tab
			
			#check if added to hits
			curquery=`echo -e "$line" | cut -f 1`
			querylength=${seclens[curquery]} #convert to length
			#querylength=`grep -m 1 "$curquery	" ./fulldb/$strainA\.fai | cut -f 2` #convert to length
			#querylength=`samtools faidx ./fulldb/$strainA $curquery | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' | cut -f 2 | tail -n 1` #convert to length
			curalnlen=`echo -e "$line" | cut -f 4`
			curquercov=`echo "$curalnlen/$querylength * 100" | bc -l | sed 's/\..*$//g'`
			if [ $curquercov -gt 39  ]; then
				hits=$((hits+1))
				echo $hits
			fi
		fi
		oldid=$curid
	done < $file
	echo "$strainA	$strainB	$hits" >> pocpnew.tmp
done


