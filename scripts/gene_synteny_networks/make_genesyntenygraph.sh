#!/bin/bash

plasmidtype="example"
touch ${plasmidtype}.network.single.out
echo -e "gene1	gene2	strain" > ${plasmidtype}.network.single.out


#for multicopy genes ie transposase just concat the locus name with clustername
for genome in `ls -1 ./*gbk`; do
	genomename=`echo $genome | sed 's/^.*\///g;s/.gbk//g'`
	
	oldlocus=""
	firstlocus="" 

	#get the number of contigs for the current input file
	numLOCUS=`grep -c "LOCUS" $genome`

	#Loop through the locus tags/genes and start outputting the syntenic gene pairs
	#Alternatively just use a file with the list of locus_tags in order instead of parsing a genbank file
	#and change the script appropriately
	for locus in `grep 'LOCUS\|/locus_tag' $genome | uniq | sed 's/^.*tag=//g;s/"//g' | sed 's/LOCUS.*/LOCUS/g'`; do
		
		#new contig, start over
		if [ "$locus" == "LOCUS" ]; then
			oldlocus=""
			continue
		fi

		#get the orthorgroup containing the current gene
		locusog=`grep "$locus	" ./ortholog_cluster_to_locus_table.txt | cut -f 1`
		locussub=`echo -e "$locus" | sed 's/_[0-9]\+$/_/g;s/_RS[0-9]\+$/_/g;s/_p[0-9]\+$/_/g'`

		#get the copy number of the current gene (does it appear more than once in an orthogroup cluster)
		ogcopynum=`grep "^$locusog" ./ortholog_cluster_to_locus_table.txt | grep -o "$locussub" | wc -l`
		
		#If the current gene is not found in any ortholog cluster use the locus tag, shouldn't happen
		if [[ -z $locusog ]]; then
			locusog=$locus
			ogcopynum=1
		fi

		#If the current gene is multi copy, output the ortholog cluster and locus tag to make a unique node
		if [ "$ogcopynum" -gt "1" ]; then
		    locusog="${locusog}_${locus}"
		fi

		#The following lines output a line describing an edge in the network. Skip the first gene but save
		#it to link to the second
		if [[ -z $firstlocus ]]; then
			firstlocus=$locusog
		fi
		if [[ ! -z $oldlocus ]]; then
			echo -e "$oldlocus	$locusog	$genomename" >> $plasmidtype\.network.single.out
		fi
		oldlocus=$locusog
	done


	#for known circular plasmids link the ends together. Remove these lines if your replicon or region is not circular
	if [[ $numLOCUS == "1" ]]; then
		echo -e "$locusog	$firstlocus	$genomename" >> $plasmidtype\.network.single.out
	fi
done
