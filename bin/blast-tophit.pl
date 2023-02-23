#!/usr/bin/perl
# ##Shawn Polson, CBCB Bioinformatics Core, University of Delaware, 2021
# FILTERS BLAST OUTPUT TO TOPHIT
# EXPECTS BLAST INPUT IN FORMAT SIMILAR TO: -outfmt "7 qaccver saccver stitle staxid sskingdoms ssciname scomname length mismatch gaps gapopen pident ppos qlen qstart qend slen sstart send frames evalue bitscore score qcovs"
# WITH COLUNS ADDING THE TAXID and TAXONOMY LINEAGE AS with taxonkit (see blast-tax.sh)



use strict;

#query acc.ver, subject acc.ver, subject title, subject tax id, subject super kingdoms, subject sci name, subject com names, alignment length, mismatches, gaps, gap opens, % identity, % positives, query length, q. start, q. end, subject length, s. start, s. end, query/sbjct frames, evalue, bit score, score, % query coverage per subject

my $infile=$ARGV[0];
my $outprefix=$infile;
$outprefix =~ s/.[^.]+$//;


#columns in btab
my $query_col="qaccver";
my $taxID_col="staxid2";
my $taxLineage_col="staxlineage";


open(my $DATA, $infile);
open(my $OUT_TOPHIT, "> ${outprefix}.tophit.btab");
open(my $OUT_VIRHIT, "> ${outprefix}.topviralhit.btab");
#open(my $OUT_ORF_SUM, "> ${outprefix}.tax.orf_summary.tsv");


my %orf_summary;
#my %results;
#my %lineages;
#my @taxids;

my $header_row=<$DATA>;
chomp $header_row;
print $OUT_TOPHIT "${header_row}\n";
print $OUT_VIRHIT "${header_row}\n";
$header_row =~ s/^[#]+//;

my @col_names=split(/\t/,$header_row);
my $curr;
my $topviral=0;

while(<$DATA>) {
	chomp;
	my $row=$_;
	my %cols;
	@cols{@col_names} = split(/\t/,$row);
	if(! $curr || $cols{$query_col} ne $curr) {
		if($curr) {
			#future additons:
			#summarize last ORF
			#LCA
			#common taxa

		}

	
		$curr=$cols{$query_col};
		print $OUT_TOPHIT "${row}\n";
		$topviral=0;

		#%lineages = ();
		#@taxids = ();
	}
	
	if ($topviral != 1 && $cols{$taxLineage_col} =~ /k__Viruses/) {
		print $OUT_VIRHIT "${row}\n";
		$topviral=1;
	}

}

close($DATA);
close($OUT_TOPHIT);
close($OUT_VIRHIT);


