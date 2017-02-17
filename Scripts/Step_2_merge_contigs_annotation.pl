#!/usr/bin/env perl

use strict;
use autodie;
require '../config.txt';
use Getopt::Long;
my $h='';
my $code_dataset='';
my $wdir='';
my $ref_phage_clusters='';
GetOptions ('help' => \$h, 'h' => \$h,'i=s'=>\$code_dataset, 'd=s'=>\$wdir, 'r=s'=>\$ref_phage_clusters);
if ($h==1 || $code_dataset eq "" || $wdir eq "" || $ref_phage_clusters eq ""){ # If asked for help or did not set up any argument
	print "# Script to generate the merged contig annotation (annotate each gene)
# Arguments :
# -i : code dataset
# -d : working directory
# -r : reference file for the viral clusters database used
";
	die "\n";
}
my $contig_list=$fastadir."/".$code_dataset."_contig_list.txt";
my %check_contig;
open my $txt,"<",$contig_list;
while(<$txt>){
	chomp($_);
	$check_contig{$_}=1;
}
close $txt;

my $circu_list=$path{'input_dir'}.$code_dataset."/".$code_dataset."_circular.txt";
my %circu;
open my $li, '<', $circu_file;
while(<$li>){
	chomp($_);
	my @tab=split("\t",$_);
	my $id_c=$tab[0];
	$circu{$id_c}=1;
}
close $li;

my $prodigal_file=$path{'input_dir'}.$code_dataset."/".$code_dataset."_Prodigal.gff";
my $n2=0;
my %size;
my %order_gene;
my %predict;
my %type;
my $id_c="";
# Read all gene predictions
open my $gff,"<",$prodigal_file;
while(<$gff>){
	chomp($_);
	if($_=~/^# Sequence Data: seqnum=\d+;seqlen=(\d+);seqhdr\"(\S+)\"/){
		if ($check_contig{$2}==1){
			$size{$2}=$1;
		}
	}
	elsif ($_=~/^#/){next}
	my @tab=split("\t",$_);
	if ($check_contig{$tab[0]}==1){
		if ($tab[8]=~/^ID=\d+_(\d+);/){
			my $id_prot=$tab[0].$1;
			$order_gene{$tab[0]}{$id_prot}=$n2;
			$predict{$tab[0]}{$id_prot}=$_;
			$n2++;
		}
		else{
			print "!!! PBLM WITH PRODIGAL LINE $_ == $tab[8]\n";
		}
	}
}
close $gff;

# first the BLAST vs unclustered , which annotation will eventually be erased by the HMM vs Phage cluster if any (that we trust more)
my $blast_vs_unclustered=$wdir.'Contigs_prots_vs_unclustered.tab';
my %affi_phage_cluster;
my $score_blast_th=50;
my $evalue_blast_th=0.001;
open my $tsv, '<', $blast_vs_unclustered;
while (<$tsv>){
	chomp($_);
	my @tab=split("\t",$_);
	my $seq=$tab[0];
	my $match=$tab[1];
	$match=~s/\|/_/g;
	my $evalue=$tab[10];
	my $score=$tab[11];
	if ($score>=$score_blast_th && $evalue<=$evalue_blast_th && (!defined($affi_phage_cluster{$seq}) || ($score>$affi_phage_cluster{$seq}{"score"})) && ($seq ne $match)){ ## We add the $seq ne $match so that we do not count a match to a phage sequence when it's only itself in the unclustered pool from a previous revision.
		$affi_phage_cluster{$seq}{"score"}=$score;
		$affi_phage_cluster{$seq}{"evalue"}=$evalue;
		$affi_phage_cluster{$seq}{"match"}=$match;
# 				print "$seq match $match\n";
	}
}
close $tsv;

my $score_th=40;
my $evalue_th=0.00001;
# Then reading the annotation from the HMM vs Viral Clusters
my $hmm_phage_clusters=$wdir.'Contigs_prots_vs_Phage_Gene_Catalog.tab';
open my $tsv, '<', $hmm_phage_clusters;
while(<$tsv>){
	chomp($_);
	if ($_=~m/^#/){
		next;
	}
	else{
		my @splign=split(m/\s+/,$_);
		my $seq=$splign[0];
		my $match=$splign[2];
		$match=~s/\.ali_faa//g;
		my $evalue=$splign[4];
		my $score=$splign[5];
		if ($score>=$score_th && $evalue<=$evalue_th && (!defined($affi_phage_cluster{$seq}) || ($score>$affi_phage_cluster{$seq}{"score"}))){
			$affi_phage_cluster{$seq}{"score"}=$score;
			$affi_phage_cluster{$seq}{"evalue"}=$evalue;
			$affi_phage_cluster{$seq}{"match"}=$match;
		}
	}
}
close $tsv;

# Then reading annotation from PFAM
my $hmm_pfama=$path{'input_dir'}.$code_dataset."/".$code_dataset."_Prodigal_Pfam.tbl";
open my $tsv, '<', $hmm_pfama;
while(<$tsv>){
	chomp($_);
	if ($_=~m/^#/){
		next;
	}
	else{
		my @splign=split(m/\s+/,$_);
		my $seq=$splign[0];
		my $match=$splign[3];
		my $evalue=$splign[6];
		my $score=$splign[7];
		if ($score>=$score_th && $evalue<=$evalue_th && (!defined($affi_pfam{$seq}) || ($score>$affi_pfam{$seq}{"score"}))){
			$affi_pfam{$seq}{"score"}=$score;
			$affi_pfam{$seq}{"evalue"}=$evalue;
			$affi_pfam{$seq}{"match"}=$match;
		}
	}
}
close $tsv;

# We also read the annotation for each phage cluster, i.e. its category
my %phage_cluster;
open my $psv, '<', $ref_phage_clusters;
while (<$psv>){
	chomp($_);
	my @tab=split(/\|/,$_);
	$phage_cluster{$tab[0]}{"category"}=$tab[1];
}
close $psv;


# Final output format:
# >Contig,nb_genes,circularity
# gene_id,start,stop,length,strand,affi_phage,score,evalue,category,affi_pfam,score,evalue,
my $out_file_affi  = catfile($wdir, $code_dataset . '_affi-contigs.csv');
open my $s1, '>', $out_file_affi;
my $n=0;
foreach(@liste_contigs){
	$n++;
	if ($n % 10000 == 0){print "$n-ieme contig\n";}
	my $contig_c=$_;
	my $circ="l";
	if ($circu{$contig_c}==1){$circ="c";}
	my @tab_genes=sort {$order_gene{$contig_c}{$a} <=> $order_gene{$contig_c}{$b} } keys %{$predict{$contig_c}};
	my $n_g=$#tab_genes+1;
	print $s1 ">$contig_c|$n_g|$circ\n";
	foreach(@tab_genes){
		my $g_c=$_;
		my @tab=split("\t",$predict{$contig_c}{$g_c});
		$g_c=$contig_c."-".$g_c;
		my $name=$tab[0];
		my $start=$tab[1];
		my $stop=$tab[2];
		my $strand=$tab[3];
		my $frame=$tab[4];
		my $affi_pc="-";
		my $affi_pc_score="-";
		my $affi_pc_evalue="-";
		my $affi_category="-";
		if (defined($affi_phage_cluster{$g_c})){
			my $phage_c=$affi_phage_cluster{$g_c}{"match"};
			if (defined($phage_cluster{$phage_c}{"category"})){$affi_category=$phage_cluster{$phage_c}{"category"};}
# 			else{print "No category for $phage_c ????????\n";} # Blast unclustered do not have any category
			$affi_pc=$phage_c;
			$affi_pc_score=$affi_phage_cluster{$g_c}{"score"};
			$affi_pc_evalue=$affi_phage_cluster{$g_c}{"evalue"};
		}
		my $affi_pfam="-";
		my $affi_pfam_score="-";
		my $affi_pfam_evalue="-";
		if (defined($affi_pfam{$g_c})){
			$affi_pfam=$affi_pfam{$g_c}{"match"};
			$affi_pfam_score=$affi_pfam{$g_c}{"score"};
			$affi_pfam_evalue=$affi_pfam{$g_c}{"evalue"};
		}
		my $length=$stop-$start;
		if ($length<0){ # It can happen if one gene overlap the contig origin
			$length=($size{$contig_c}-$start)+$stop;
		}
		print $s1 "$g_c|$start|$stop|$length|$strand|$affi_pc|$affi_pc_score|$affi_pc_evalue|$affi_category|$affi_pfam|$affi_pfam_score|$affi_pfam_evalue\n";
	}
}
close $s1;
