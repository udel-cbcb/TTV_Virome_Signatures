#!/usr/bin/perl
##Shawn Polson, CBCB Bioinformatics Core, University of Delaware, 2021
# OUTPUTS CONSENSUS CONTIG TAXONOMY INFORMATION
# EXPECTS B:AST INPUT IN FORMAT SIMILAR TO: -outfmt "7 qaccver saccver stitle staxid sskingdoms ssciname scomname length mismatch gaps gapopen pident ppos qlen qstart qend slen sstart send frames evalue bitscore score qcovs"
# WITH COLUNS ADDING THE TAXID and TAXONOMY LINEAGE AS with taxonkit (see blast-tax.sh)

#query acc.ver, subject acc.ver, subject title, subject tax id, subject super kingdoms, subject sci name, subject com names, alignment length, mismatches, gaps, gap opens, % identity, % positives, query length, q. start, q. end, subject length, s. start, s. end, query/sbjct frames, evalue, bit score, score, % query coverage per subject

use strict;

my $infile=$ARGV[0];
my $outprefix=$infile;
$outprefix =~ s/.[^.]+$//;
our $maj_pct_cutoff=0.50;
my $maj_pct_co_format=$maj_pct_cutoff;
$maj_pct_co_format=~s/^[01].//;

#columns in btab
our $query_col="qaccver";
our $subject_col="saccver";
our $contig_col="qcontig";
our $taxid_col="staxid2";
our $lineage_col="staxlineage";


open(my $DATA, $infile);
open(our $OUT_CTG_SUM, "> ${outprefix}.tax_ctg_summary.tsv");
print $OUT_CTG_SUM "contig\tlca_taxid\tlca_lineage\tcommon_taxid\tcommon_lineage\tcommon_pct\tgte${maj_pct_co_format}_taxid\tgte${maj_pct_co_format}_lineage\tgte${maj_pct_co_format}_pct\n";


my @lineages;
my @taxids;

my $header_row=<$DATA>;
chomp $header_row;
$header_row =~ s/^[#]+//;

my @col_names=split(/\t/,$header_row);
my $curr;

while(<$DATA>) {
	chomp;
	my $row=$_;
	my %cols;
	@cols{@col_names} = split(/\t/,$row);
	if(! $curr || $cols{$contig_col} ne $curr) {
		if($curr) {
			process_contig($curr, \@taxids, \@lineages);
		}

	
		$curr=$cols{$contig_col};

		@lineages = ();
		@taxids = ();
	}
	push(@lineages, $cols{$lineage_col});
	push(@taxids, $cols{$taxid_col});
}
###NEED TO PROCESS LAST CONTIG
process_contig(\@taxids,\@lineages);


close($DATA);
close($OUT_CTG_SUM);

sub process_contig {
	my $curr=shift;
	my @taxids=@{ shift @_ };
	my @lineages=@{ shift @_ };
	
	#LCA
	my $taxid_list=join(" ",@taxids);
	#print $taxid_list."\n";
	my @lca=split(/\t/,`echo "$taxid_list" | taxonkit lca --threads 4| taxonkit reformat --taxid-field 2 --format "{k};{p};{c};{o};{f};{g};{s};{t}" --fill-miss-rank --pseudo-strain --miss-rank-repl-prefix "norank_" --trim --output-ambiguous-result --add-prefix --threads 4`);
	my $lca_taxid=$lca[1];
	my $lca_lineage=$lca[2];
	chomp $lca_lineage;

	#most common lineage
	my %taxid_hash;
	foreach my $i (@taxids) {
        	$taxid_hash{$i}++;
	}
	my $common_taxid=(sort {$taxid_hash{$b} <=> $taxid_hash{$a}} keys %taxid_hash)[0];
	my $common_pct=$taxid_hash{$common_taxid}/($#taxids+1); 
	my $common_lineage=(split(/\t/,`echo "$common_taxid" | taxonkit reformat --taxid-field 1 --format "{k};{p};{c};{o};{f};{g};{s};{t}" --fill-miss-rank --pseudo-strain --miss-rank-repl-prefix "norank_" --trim --output-ambiguous-result --add-prefix --threads 4`))[1];
	chomp $common_lineage;

	#majority level
	my $maj_taxid="";
	my $maj_pct="";
	my $maj_lineage="";
	if($common_pct >= $maj_pct_cutoff) {
		$maj_taxid=$common_taxid;
		$maj_pct=$common_pct;
		$maj_lineage=$common_lineage;
	} else {
		foreach (1..7) {
			my @tmp_lineages;
			foreach my $j (@lineages) {
				$j =~ s/;[^;]+$//;
				push(@tmp_lineages, $j);
			}

			my %tmp_lineage_hash;
			foreach my $j (@tmp_lineages) {
				$tmp_lineage_hash{$j}++;
			}
			my $tmp_maj_lineage=(sort {$tmp_lineage_hash{$b} <=> $tmp_lineage_hash{$a}} keys %tmp_lineage_hash)[0];
			my $tmp_maj_pct=$tmp_lineage_hash{$tmp_maj_lineage}/($#tmp_lineages+1);
			#print $tmp_maj_lineage."\t".$tmp_maj_pct."\n";
			if($tmp_maj_pct >= $maj_pct_cutoff) {
				$maj_lineage=$tmp_maj_lineage;
				chomp $maj_lineage;
				$maj_pct=$tmp_maj_pct;
				$tmp_maj_lineage =~ s/[kpcofgst]__//g;
				$tmp_maj_lineage =~ s/;$//;
				#print $tmp_maj_lineage."\n";
				my $maj_taxid_string=(split(/\t/,`echo "$tmp_maj_lineage" | taxonkit reformat --threads 4 --show-lineage-taxids --lineage-field 1`))[-1];
				$maj_taxid=(split(/;/, $maj_taxid_string))[-1];
				chomp $maj_taxid;
				last;
			}
		}
	}

        
	#summarize last contig
	print $OUT_CTG_SUM "${curr}\t${lca_taxid}\t${lca_lineage}\t${common_taxid}\t${common_lineage}\t${common_pct}\t${maj_taxid}\t${maj_lineage}\t${maj_pct}\n";

}
