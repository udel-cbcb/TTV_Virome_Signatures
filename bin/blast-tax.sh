#!/bin/bash

#Shawn Polson, CBCB Bioinformatics Core, University of Delaware, 2021
#REQUIRES taxonkit and csvtk be in system path
#INPUT: 1 or more blast results from blastp or RUBBLE formatted:
#   -outfmt "7 qaccver saccver stitle staxid sskingdoms ssciname scomname length mismatch gaps gapopen pident ppos qlen qstart qend slen sstart send frames evalue bitscore score qcovs"


THREADS=       # NUMBER OF THREADS TO PARALLELIZE
INFOFN=        # Path to a File of file paths to analyze (1 per line) with blast results in 
MAPDIR=        # BBMAP Reference directory
OUTDIR=        # Output directory
OUTFOFN=       # Output file of file paths

# COMMAND(s) TO RUN
mkdir -p ${OUTDIR}
cd ${OUTDIR}

for i in $(cat ${INFOFN})
do
	PREFIX=$(basename ${i} .btab)
	MAPPREFIX=$(echo $PREFIX|sed -E -e "s/^([^.]+)\..+$/\1/")
	TMPFILE=${OUTDIR}/${PREFIX}.tmptax.btab
	OUTFILE=${OUTDIR}/${PREFIX}.wtax.btab
	egrep -v "^#" $i | sed -E -e "s|TaxID=N/A|TaxID=1|" -e "s/^(.+TaxID=([0-9]+) .+)$/\1	\2/" > ${TMPFILE}
	taxonkit reformat --taxid-field 26 --format "{k};{p};{c};{o};{f};{g};{s};{t}" --fill-miss-rank --pseudo-strain --miss-rank-repl-prefix "norank_" --trim --output-ambiguous-result --add-prefix --threads ${THREADS} ${TMPFILE} > ${OUTFILE}
	rm ${TMPFILE}
	sed -Ei "1i $(head -n 1 $i | sed -E -e "s/$/	staxid2	staxlineage/")" ${OUTFILE}
	blast-tophit.pl ${OUTFILE}
	blast-tax-contig.pl ${OUTDIR}/${PREFIX}.wtax.tophit.btab
	csvtk join -t -f 1 --left-join -C % ${OUTDIR}/${PREFIX}.wtax.tophit.tax_ctg_summary.tsv ${MAPDIR}/${MAPPREFIX}.constats.txt > ${OUTDIR}/${PREFIX}.wtax.tophit.tax_ctg_summary.wcov.tsv
	blast-tax-contig.pl ${OUTDIR}/${PREFIX}.wtax.topviralhit.btab
	csvtk join -t -f 1 --left-join -C % ${OUTDIR}/${PREFIX}.wtax.topviralhit.tax_ctg_summary.tsv ${MAPDIR}/${MAPPREFIX}.constats.txt > ${OUTDIR}/${PREFIX}.wtax.topviralhit.tax_ctg_summary.wcov.tsv	
	realpath ${OUTFILE} >> ${OUTFOFN}
done

wait
