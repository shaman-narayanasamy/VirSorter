#!/usr/bin/env perl

use strict;
use autodie;
require '../config.txt';
use Getopt::Long;
my $h='';
my $code_dataset='';
my $fastadir='';
my $th_nb_genes='';
GetOptions ('help' => \$h, 'h' => \$h,'i=s'=>\$code_dataset,'f=s'=>\$fastadir,'t=s'=>\$th_nb_genes);
if ($h==1 || $code_dataset eq "" || $fastadir eq "" || $th_nb_genes eq ""){ # If asked for help or did not set up any argument
	print "# Script to select relevant contigs (based on # of predicted ORFs) and get a nice fasta file of predicted proteins
# NOTE THAT WE USE PRODIGAL PREDICTIONS - NO PARTIAL AT ALL (?)
# Arguments :
# -i : code dataset
# -f : path to fasta dir for VirSorter
# -t : threshold on gene numbers
";
	die "\n";
}

my $prodigal_file=$path{'input_dir'}.$code_dataset."/".$code_dataset."_Prodigal.gff";
my %count;
open my $gff,"<",$prodigal_file;
while(<$gff>){
	chomp($_);
	if ($_=~/^#/){next;}
	my @tab=split("\t",$_);
	$count{$tab[0]}++;
}
close $gff;

my $faa_file=$path{'input_dir'}.$code_dataset."/".$code_dataset."_Prodigal_proteins.faa";
my $prot_file=$fastadir."/".$code_dataset."_prots.fasta";
my $contig_list=$fastadir."/".$code_dataset."_contig_list.txt";
open my $faa,"<",$faa_file;
open my $s1,">",$prot_file;
open my $s2,">",$contig_list;
my $tag=0;
my %seen;
while(<$faa>){
	chomp($_);
	if ($_=~/^>(\S+)/){
		my $id=$1;
		$tag=0;
		if ($id=~/(.*)_\d+/){
			my $id_contig=$1;
			if ($count{$id_contig}>=$th_nb_genes){
				if ($seen{$id_contig}==1){}
				else{
					print $s2 "$id_contig\n";
					$seen{$id_contig}=1;
				}
				print $s1 ">$id\n";
				$tag=1;
			}
		}
		else{
			print "!! PBLM WITH ID $id\n";
		}
	}
	elsif($tag==1){
		print $s1 "$_\n";
	}
}
close $faa;
close $s1;
close $s2;