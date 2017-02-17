#!/usr/bin/env perl
=head1 SYNOPSIS
	wrapper_phage_contigs_sorter_iPlant.pl -s spid
Required Arguments:
	-s	spid
	--help	Show help and exit

=head1 DESCRIPTION

Wrapper for detection of viral contigs

=cut

## NEED TO LOAD THE FOLLOWING MODULES ###
## module load blast
## module load hmmer

use strict;
use warnings;
require 'config.txt';
my $help='';
my $n_cpus=16;
my $data_dir=$path{'data_dir'};
my $wdir = $path{'w_dir'};
my $Bin = $path{'bin_dir'};
my $code_dataset='';
GetOptions('s=s'=> \$code_dataset);
my $input_file=$path{'input_dir'}.$code_dataset."/".$code_dataset."_contigs.fasta"; ## NOTE THIS COULD CHANGE TO .FNA ##
my $choice_database=2;
if ($help) {
    pod2usage();
}
unless ($input_file) {
    pod2usage('Missing FASTA file');
}

say map { sprintf "%-15s: %s\n", @$_ } (
	['Bin',		$Bin],
	['Dataset',	$code_dataset],
	['Input file',	$input_file],
	['Db',		$choice_database],
	['Working dir',	$wdir],
);

my $path_hmmsearch     = which('hmmsearch') or die "Missing hmmsearch\n";
my $path_last        = which('blastp')    or die "Missing blastp\n";
my $script_dir         = catdir($Bin, 'Scripts');
my $readme_file        = catfile($data_dir, 'VirSorter_Readme.txt');

my $generic_ref_file = catfile($data_dir,'Generic_ref_file.refs');

my $dir_Phage_genes    = catdir($data_dir,'Db/Phage_gene_catalog');
my $ref_phage_clusters = catfile($data_dir, 'Db/Phage_gene_catalog', 'Phage_Clusters_current.tab');
if ($choice_database == 2) {
    $dir_Phage_genes    = catdir($data_dir, 'Db/Phage_gene_catalog_plus_viromes');
    $ref_phage_clusters = catfile($data_dir,'Db/Phage_gene_catalog_plus_viromes', 'Phage_Clusters_current.tab');
}
my $out = "";
## SETTING UP THE WORKING DIRECTORY
my $log_dir = catdir($wdir, 'logs');
mkpath($log_dir);
my $log_out = catfile($log_dir, 'out');
my $log_err = catfile($log_dir, 'err');
die("");
# cp fasta file in the wdir
my $fastadir = catdir($wdir, 'fasta');
if ( !-d $fastadir ) {
    mkpath($fastadir);
    # NEW STEP 1 -> SELECT ONLY CONTIGS OF 2 ORFS OR MORE, AND EXTRACT THEIR PROTEIN
    my $nb_gene_th = 2; # At least two complete genes on the contig
    &run_cmd($script_dir."Step_1_contigs_cleaning_and_gene_prediction.pl -i $code_dataset -f $fastadir -t $nb_gene_th >> $log_out 2>> $log_err");
}

my $fasta_file_prots = catfile($fastadir, $code_dataset . "_prots.fasta");

# # Files that will stay along the computations
my $out_hmmsearch = catfile($wdir, 'Contigs_prots_vs_Phage_Gene_Catalog.tab');
my $out_hmmsearch_bis = catfile($wdir, 'Contigs_prots_vs_Phage_Gene_Catalog.out');
my $out_blast_unclustered = catfile($wdir, 'Contigs_prots_vs_Phage_Gene_unclustered.tab');
my $out_file_affi  = catfile($wdir, $code_dataset . '_affi-contigs.csv');
my $out_file_phage_fragments = catfile($wdir, $code_dataset . '_phage-signal.csv');
my $global_out_file = catfile($wdir, $code_dataset . '_global-phage-signal.csv');
my $new_prots_to_cluster = catfile($wdir, $code_dataset . '_new_prot_list.csv');
# # Get the final result file ready
&run_cmd("touch $global_out_file");
my $r_n = -1;
# Si on a des nouvelles prots a clusteriser ou si on est dans la premiere
# revision
#
while ( (-e $new_prots_to_cluster || $r_n == -1) && ($r_n<=10) ) {
	$r_n++;    # New revision of the prediction
	my $dir_revision = catdir($wdir, 'r_' . $r_n);
	print "### Revision $r_n\n";
	if (!-d $dir_revision) {
		## mkdir for this revision
		mkpath($dir_revision);
		if ( $r_n == 0 ) {
			# First revision, we just import the Refseq database
			mkpath(catdir($dir_revision, 'db'));
			### WE COMMENT THIS, BECAUSE WE SHOULD NEVER HAVE CUSTOM PHAGE ADDITION IN AUTOMATIC PROCESSING, BUT LEAVE IT THERE JUST IN CASE
#			## Adding custom sequences to the database if required by the user
#			if ($custom_phage ne '') {
#				my $script_custom_phage = catfile(
# 				$script_dir, "Step_first_add_custom_phage_sequence.pl"
#				);
# 				$out = `$script_custom_phage $custom_phage $dir_Phage_genes/ $dir_revision/db >> $log_out 2>> $log_err`;
#				say "Adding custom phage to the database : $out";
#			}
# 			# should replace Pool_cluster / Pool_unclustered and
#			# Pool_new_unclustered else , we just import the Refseq database
#			else { 
			&run_cmd("cp $dir_Phage_genes/* $dir_revision/db/"); 
# 			}
		}
		else {
			my $previous_r = $r_n - 1;
			my $previous_fasta_unclustered = catfile($wdir, 'r_'. $previous_r, 'db', 'Pool_unclustered.faa');
			print "Step 1.1 new clusters and new database";
			&run_cmd("$script_dir/Step_0_make_new_clusters.pl $dir_revision $fasta_file_prots $previous_fasta_unclustered $new_prots_to_cluster >> $log_out 2>> $log_err");
			# Rm the list of prots to be clustered now that they should be clustered
			&run_cmd("rm $new_prots_to_cluster");
		}
		# Check if there are some data in these new clusters, or if all the new proteins are unclustered
		my $new_db_profil = catfile($dir_revision, 'db', 'Pool_clusters.hmm');
		my $check = 0;
		if (-s $new_db_profil) {
			open my $DB, '<', $new_db_profil;
			while (<$DB>) {
				chomp($_);
				if ( $_ =~ /^NAME/ ) { 
					$check++; 
				}
			}
			close $DB;
		}
		if ($check == 0) {
			print "There are no clusters in the database, so skip the hmmsearch\n";
		}
		else {
			my $out_hmmsearch_new =catfile($dir_revision, 'Contigs_prots_vs_New_clusters.tab');
			my $out_hmmsearch_bis_new =catfile($dir_revision, 'Contigs_prots_vs_New_clusters.out');
			print "Step 1.2 : HMMSEARCH against the viral database \n";
			&run_cmd("$path_hmmsearch --tblout $out_hmmsearch_new --cpu $n_cpus -o $out_hmmsearch_bis_new --noali $new_db_profil $fasta_file_prots >> $log_out 2>> $log_err");
			&run_cmd("cat $out_hmmsearch_new >> $out_hmmsearch");
		}
		my $out_blast_new_unclustered =catfile($dir_revision, 'Contigs_prots_vs_New_unclustered.tab');
		my $blastable_unclustered =catfile( $dir_revision, 'db', 'Pool_new_unclustered' );
		print "Step 1.3 :  BLAST against viral gene singletons (i.e. not in HMM profiles)";
		&run_cmd("$path_blastp -query $fasta_file_prots -db $blastable_unclustered -out $out_blast_new_unclustered -num_threads $n_cpus -outfmt 6 -evalue 0.001 >> $log_out 2>> $log_err");
		&run_cmd("cat $out_blast_new_unclustered >> $out_blast_unclustered");
		# Make backup of the previous files to have 
		# trace of the different steps
		my $backup_affi = catfile($dir_revision, 'affi_backup.csv');
		my $backup_phage_signal = catfile($dir_revision, 'phage_signal_backup.csv');
		my $backup_global_signal = catfile($dir_revision, 'global_signal_backup.csv');
		if (-e $out_file_affi) {`cp $out_file_affi $backup_affi`;}
		if (-e $out_file_phage_fragments) {`cp $out_file_phage_fragments $backup_phage_signal`;}
		if (-e $global_out_file) {`cp $global_out_file $backup_global_signal`;}
	}
	## Complete the affi
	print "Step 2 : Merging annotations \n";
	&run_cmd("$script_dir/Step_2_merge_contigs_annotation.pl -i $code_dataset -d $wdir -r $ref_phage_clusters >> $log_out 2>> $log_err");
	## This generate a csv table including the map of each contig, with PFAM and Viral PCs annotations, as well as strand and length of genes
	print "Step 3 : Detecting viral contigs and/or regions\n";
	&run_cmd("$script_dir/Step_3_highlight_phage_signal.pl -i $code_dataset -d $wdir -r $generic_ref_file >> $log_out 2>> $log_err");
	# Decide which contigs are entirely viral and which are prophages, and
	# which of both of these categories are phage enough to be added to the
	# databases
	print "Step 4: Setting up the final result file - $global_out_file\n";
	&run_cmd("$script_dir/Step_4_summarize_phage_signal.pl -i $code_dataset -d $wdir -n $new_prots_to_cluster >> $log_out 2>> $log_err");
}

print "Step 5 : Prepping output files\n";
my $final_affi_file=catfile($path{'input_dir'}, $code_dataset . '_Prodigal_VirSorter_Annotation.psv');
&run_cmd("cp $out_file_affi $final_affi_file");

my $final_result_file=catfile($path{'input_dir'}, $code_dataset . '_Prodigal_VirSorter_Result.csv');
&run_cmd("cp $out_file_phage_fragments $final_result_file");

my $final_result_file=catfile($path{'input_dir'}, $code_dataset . '_Prodigal_VirSorter_viruses.tab');
open my $s1,">",$final_result_file;
print $s1 "# Contig\tConfidence (1 sure, 2 pretty sure, 3 not so sure)\tType (all viral contig or prophage)\tCoordinates (if prophage)\n";
open my $csv,"<",$out_file_phage_fragments;
my $c_i=0;
while(<$csv>){
	chomp($_);
	if ($_=~/^## (\d)/){$c_i=$1;}
	elsif($_=~/^##/){next}
	my @tab=split(",",$_);
	if ($c_i<=3){
		print $s1 "$tab[0]\t$c_i\tall viral\t\n";
	}
	else{
		my $cat=$c_i-3;
		$tab[2]=~/.*-(gene_\d+-gene_\d+)/;
		print $s1 "$tab[0]\t$c_i\tprophage\t$1\n";
	}
}
close $csv;
close $s1;

my $datestring = localtime();
my $virsorter_readme=catfile($path{'input_dir'}, $code_dataset . '_Prodigal_VirSorter_readme.txt');
open my $s1, '>', $virsorter_readme;
print $s1 "VirSorter parameters used :\n";
print $s1 "--> Fasta file mined for viral sequences : $input_file\n";
print $s1 "--> Viral database used : \n";
if ($choice_database == 2) {
    say $s1 join(' ', 
        "Viromes : all bacterial and archaeal virus genomes in Refseq,",
        "as of January 2014, plus non-redundant predicted genes from viral",
        "metagenomes (including seawater, freshwater, and human-related",
        "samples)"
    );
}
else {
    print $s1 "RefseqABVir (all bacterial and archaeal virus genomes " .
        "in Refseq, as of January 2014)\n";
}
print $s1 join(' ',
        "VirSorter was run with the in the 'Virome Decontamination' mode:",
        "overall metrics for microbial sequences were not evaluated from the",
        "complete dataset, but instead pre-computed values based on bacterial",
        "and archaeal genomes from Refseq were used."
    )."\n";

print $s1 "This VirSorter computation finished on $datestring\n";
close $s1;

# Plus clean the output directory
print "Step 6: Cleaning the output directory\n";
# &run_cmd("rm -r $wdir");



sub run_cmd{
	my $cmd=$_[0];
	print "$cmd\n";
	my $out=`$cmd`;
	print "$out\n";
}