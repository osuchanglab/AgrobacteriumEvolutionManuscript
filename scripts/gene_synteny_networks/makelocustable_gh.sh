#!/bin/bash

#Example code to create a locus tag/ortholog cluster table from get_homologues
#cluster files. Change "mydataset" to the name of your folder containing the 
#fasta file cluster output. 
for file in `ls -1 ./mydataset_homologues/mydataset_prot_OMCL/*faa`; do
	dataset=`echo $file | sed 's/^.*\///g;s/\.faa//g'`
	loci=`grep '>' $file | sed 's/ |.*//g;s/>ID://g' | tr '\n' '\t'`
	echo "$dataset	$loci"
done


#Alternatively transpose the "pangenome_matrix_genes_t0.tab" output of get_homologues 
awk '
{
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' pangenome_matrix_genes_t0.tab > pangenome_matrix_genes_t0.tr.tab

