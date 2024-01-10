#!/use/bin/perl
use strict;
use warnings;

# Example usage: perl produce_generic_exomiser_master_file_final.pl exomiser_samples_b38.tsv denovo_b38.tsv ccr_95_hg38.tsv b38

my $samples_file = $ARGV[0] || die;# list of sample IDs used to name Exomiser output files (one per line)
my $de_novo_file = $ARGV[1] || '';# optional tab-sep file with sampleID, chr, pos, ref, alt for de novo called variants 
my $ccr_file = $ARGV[2] || '';# optional CCR file defining constrained coding regions in chr:start-end format 
my $assembly = $ARGV[3] || 'b38';# b37 or b38(default) assembly

my %de_novo_new_calls;
if ($de_novo_file ne ''){
    open(IN1,"$de_novo_file");
    while (my $line = <IN1>){
	chop $line;
	my @line = split(/\t/,$line);
	my $sample = $line[0];
	my $key = $line[1].":".$line[2].":".$line[3].":".$line[4];
	$de_novo_new_calls{$sample}{$key} = 'Y';
    }
}
close IN1;

my %ccr;
if ($ccr_file ne ''){
    open(IN,"$ccr_file");
    while (my $line= <IN>){
	chop $line;
	my ($chr,$range) = split(/\:/,$line);
	$chr =~ s/chr//g;
	my ($start,$end) = split(/\-/,$range);
	for(my $pos = $start; $pos<= $end; $pos++){
	    $ccr{$chr}{$pos} = 1;
	}
    }
    close IN;
}

my $total_case_count = 0;
my %var_counts;

# add fam size, hpo ids and hpo terms to the samples file and add to the output at the end

open(EXOMISER_MASTER_FILE,">exomiser_master_file_final_allvars.tsv");
print EXOMISER_MASTER_FILE "Sample ID\tGene\tVariant\tGenotypes\tExomiserMOI\tFunctional Class\tMax Freq\tGnomAD Freq\tDe novo\tRank\tScore\tVariant Score\tPheno Score\tHuman Pheno Score\tBest Evidence\tHuman Evidence\tMouse Evidence\tFish Evidence\tHuman PPI Evidence\tHGVS\tAssembly\tFam Structure\tHPO IDs\tHPO terms\n";
open(IN,"$samples_file");
while (my $sample_line = <IN>){
    chop $sample_line;
    my ($exomiser_sample,$fam_structure,$hpo_ids,$hpo_terms) = split(/\t/,$sample_line);
    $total_case_count++;
    my %exomiser_results;
    my %func_class_store;
    my %genotype_store;
    my %hgvs_store;
    my %exomiser_results_moi;
    my %pheno_scores;
    my %human_pheno_scores;
    my %best_evidence;
    my %human_evidence;
    my %mouse_evidence;
    my %fish_evidence;
    my %human_ppi_evidence;
    my %max_freq_store;
    my %gnomad_freq;

    my $exomiser_gene_file = $exomiser_sample.".genes.tsv";
    open(IN2,"$exomiser_gene_file");
    my $line = <IN2>;
    while ($line = <IN2>){
        chop $line;
        my @line = split(/\t/,$line);
        my $moi = $line[4];
        my $gene = $line[2];
        my $pheno_score = $line[7];
        my $human_score = $line[9];
        my $mouse_score = $line[10];
        my $fish_score = $line[11];
        my $ppi_score = $line[12];
        if (!$pheno_scores{$gene} || $pheno_score > $pheno_scores{$gene}){
            $pheno_scores{$gene} = $pheno_score;
            $human_pheno_scores{$gene} = $human_score;
            $human_evidence{$gene} = $line[16];
            $mouse_evidence{$gene} = $line[17];
            $fish_evidence{$gene} = $line[18];
            $human_ppi_evidence{$gene} = $line[19];
            if ($line[13] == $human_score){
                $best_evidence{$gene} = $line[16];
            }
            elsif ($line[13] == $mouse_score){
                $best_evidence{$gene} = $line[17];
            }
            elsif ($line[13] == $fish_score){
                $best_evidence{$gene} = $line[18];
            }
            elsif ($line[13] == $ppi_score){
                $best_evidence{$gene} = $line[19] || $line[20] || $line[21];
            }
        }
    }
    close IN2;
    
    my $exomiser_file = $exomiser_sample.".variants.tsv";
    open(IN2,"$exomiser_file");
    $line = <IN2>;
    my $last_gene = '';
    my $last_gene_combined_score = '';
    my $rank = '';
    while ($line = <IN2>){
        chop $line;
        my @line = split(/\t/,$line);
	my $moi = $line[4];
        my $chromosome = $line[14];
        my $position = $line[15];
        my $reference = $line[17];
        my $alternate = $line[18];
        my $gene = $line[2];
        my $gene_variant_score = $line[8];
	my $key = "$chromosome:$position:$reference:$alternate";
	$func_class_store{$key} = $line[23];
	$hgvs_store{$key} = $line[24];
	$genotype_store{$key} = $line[22];
	$max_freq_store{$key} = $line[36];
        my $all_freq = $line[37];
        $gnomad_freq{$key} = $all_freq;

	if ($last_gene ne '' && $gene ne $last_gene){
	    my ($de_novo,$gene_var_score,$exomiser_results_moi,$func_class,$max_freq,$genotypes,$gnomad_data,$hgvs_data,$path_data) = ('','','','','','','','','');
	    foreach my $key (sort keys %{$exomiser_results{$last_gene}}){
		$var_counts{$key}{$exomiser_sample} = 1;
		if ($de_novo_new_calls{$exomiser_sample}{$key}){
		    $de_novo = $de_novo_new_calls{$exomiser_sample}{$key}."|".$de_novo;
		}
		$func_class = $func_class_store{$key}."|".$func_class;
		$genotypes = $genotype_store{$key}."|".$genotypes;
		$hgvs_data = $hgvs_store{$key}."|".$hgvs_data;
		$exomiser_results_moi = $exomiser_results_moi{$last_gene}{$key};
                $max_freq = $max_freq_store{$key}."|".$max_freq;
                $gnomad_data =  $gnomad_freq{$key};
		$exomiser_results_moi = $exomiser_results_moi{$last_gene}{$key};
		$gene_var_score = $exomiser_results{$last_gene}{$key};
	    }
	    my $cand_vars = join("|",sort keys %{$exomiser_results{$last_gene}});
	    my $best_evidence = $best_evidence{$last_gene} || '';
            my $human_evidence = $human_evidence{$last_gene} || '';
            my $mouse_evidence = $mouse_evidence{$last_gene} || '';
            my $fish_evidence = $fish_evidence{$last_gene} || '';
            my $human_ppi_evidence = $human_ppi_evidence{$last_gene} || '';
	    print EXOMISER_MASTER_FILE "$exomiser_sample\t$last_gene\t$cand_vars\t$genotypes\t$exomiser_results_moi\t$func_class\t$max_freq\t$gnomad_data\t$de_novo\t$rank\t$last_gene_combined_score\t$gene_var_score\t".$pheno_scores{$last_gene}."\t".$human_pheno_scores{$last_gene}."\t".$best_evidence."\t$human_evidence\t$mouse_evidence\t$fish_evidence\t$human_ppi_evidence\t$hgvs_data\t$assembly\t$fam_structure\t$hpo_ids\t$hpo_terms\n";
	}

	$last_gene = $gene;		
        $exomiser_results{$gene}{$key} = $gene_variant_score if ($line[10] == 1);
	$last_gene_combined_score = $line[6];	
	$exomiser_results_moi{$gene}{$key} = $moi if ($line[10] == 1);
	$rank = $line[0] if ($line[10] == 1);
    }
    # print out the last gene result from the file that would be missed otherwise - repeat of the code above so could be a subroutine
    my ($de_novo,$gene_var_score,$exomiser_results_moi,$func_class,$max_freq,$genotypes,$gnomad_data,$hgvs_data,$path_data) = ('','','','','','','','','');
    foreach my $key (sort keys %{$exomiser_results{$last_gene}}){
	if ($de_novo_new_calls{$exomiser_sample}{$key}){
	    $de_novo = $de_novo_new_calls{$exomiser_sample}{$key}."|".$de_novo;
	}
	$func_class = $func_class_store{$key}."|".$func_class;
	$genotypes = $genotype_store{$key}."|".$genotypes;
	$hgvs_data = $hgvs_store{$key}."|".$hgvs_data;
	$exomiser_results_moi = $exomiser_results_moi{$last_gene}{$key};
	$max_freq = $max_freq_store{$key}."|".$max_freq;
	$gnomad_data =  $gnomad_freq{$key};
	$gene_var_score = $exomiser_results{$last_gene}{$key};
    }
    my $cand_vars = join("|",sort keys %{$exomiser_results{$last_gene}});
    my $best_evidence = $best_evidence{$last_gene} || '';
    my $human_evidence = $human_evidence{$last_gene} || '';
    my $mouse_evidence = $mouse_evidence{$last_gene} || '';
    my $fish_evidence = $fish_evidence{$last_gene} || '';
    my $human_ppi_evidence = $human_ppi_evidence{$last_gene} || '';
    print EXOMISER_MASTER_FILE "$exomiser_sample\t$last_gene\t$cand_vars\t$genotypes\t$exomiser_results_moi\t$func_class\t$max_freq\t$gnomad_data\t$de_novo\t$rank\t$last_gene_combined_score\t$gene_var_score\t".$pheno_scores{$last_gene}."\t".$human_pheno_scores{$last_gene}."\t".$best_evidence."\t$human_evidence\t$mouse_evidence\t$fish_evidence\t$human_ppi_evidence\t$hgvs_data\t$assembly\n";
    
    close IN2;
}

close IN;
close EXOMISER_MASTER_FILE;
# removed non-PASS variants
open(EXOMISER_MASTER_FILE_1,"exomiser_master_file_final_allvars.tsv");
open(EXOMISER_MASTER_FILE_FINAL,">exomiser_master_file_final_passvars.tsv");
my $line = <EXOMISER_MASTER_FILE_1>;
chop $line;
print EXOMISER_MASTER_FILE_FINAL "$line\tExomiser result count\tFreq in Exomiser result\tCCR Flag\n";
while ($line = <EXOMISER_MASTER_FILE_1>){
    chop $line;
    my @line = split(/\t/,$line);
    my $var = $line[2];
    my $case_count_summary = '';
    my $var_freq_summary = '';
    my $pass = 'PASS';
    my $in_ccr = 0;
    my @var = split(/\|/,$var);
    foreach my $v(@var){
        my $case_count = scalar(keys %{$var_counts{$v}});
	my $var_freq = 100*$case_count/$total_case_count;
        $case_count_summary = $case_count_summary.$case_count."|";
        $var_freq_summary = $var_freq_summary.$var_freq."|";
        if (scalar(@var) > 1){
            $pass = 'FAIL' if $var_freq > 2;
        }
        elsif ($v =~ /M/){
            $pass = 'FAIL' if $var_freq > 0.2;
        }
        else{
            $pass ='FAIL' if $var_freq > 0.1;
        }
	my ($chr,$pos,$ref,$alt) = split(/\:/,$var);
        if ($ccr{$chr}{$pos}){
            $in_ccr = 1;
        }
    }
    print EXOMISER_MASTER_FILE_FINAL "$line\t$case_count_summary\t$var_freq_summary\t$in_ccr\n" if $pass eq 'PASS';

}

close EXOMISER_MASTER_FILE_1;
close EXOMISER_MASTER_FILE_FINAL;
