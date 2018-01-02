DEVELOPMENT_PROGRAMS := \
	add_abspos \
	analyze_transmission \
	anchor_reads \
	anchor_repeatness \
	annotate_repeats \
	assess_known_cn_events \
	bam_genotyper \
	bin_events \
	bin_p \
	bridge_counts \
	bridge_effect \
	bridge_gene2 \
	bridge_gene \
	cancer_pop \
	candidate_stats \
	check_bridges \
	check_bridges_result \
	check_cand_cn \
	check_for_pseudogene \
	check_mappability \
	check_mumdex \
	combine_finebins \
	compare_segmentation \
	compare_yoonha \
	control_mview \
	count_anchors \
	count_genes \
	count_pseudogenes \
	debruijn \
	denovo_cn \
	denovo_pseudogenes \
	determine_cn_sex \
	determine_sex \
	dsDNAvt \
	encode \
	event_histogram \
	family_count \
	family_figure \
	fastq_kmer_counter \
	fastq_mapper \
	fastq_vt_mapper \
	find_bridge \
	find_indel \
	find_microsatellite \
	gene_view \
	jackpot \
	karyotype \
	local_similarity \
	lookup_bridge \
	male_cn_stats \
	mapper \
	matrix \
	missing_chromosome \
	model_transmission \
	mumdex_examples \
	mumdex_sequences \
	museq_primers \
	optional_test \
	pair_view \
	pivot_table \
	plot_cn_detailed \
	plot_cn_stats \
	plot_coverage \
	plot \
	plot_invariants \
	plot_simulate_cn \
	plot_tnp \
	plot_unequal \
	poisson \
	position_coverage \
	primers \
	print_invariants \
	process_bisulfite_rna \
	random_mumdex_sequences \
	random_position \
	rare \
	repeatness \
	repeats \
	show_all_counts \
	show_counts \
	show_counts_special \
	show_mums \
	show_pairs \
	similar_bridge \
	simulate_cn \
	simulate_events2 \
	simulate_events \
	skbr3_in_skn1 \
	smooth_data \
	snp_report \
	sssa \
	subsample \
	test_mumdex \
	test_numerical \
	test_psplot \
	test_threads \
	test_x11 \
	transmission \
	unequal_bridges \
	validate_leukemia \
	validate_mumdex \
	wg_smash_comparison \
	x11plot

PROGRAMS += $(DEVELOPMENT_PROGRAMS)
all : $(DEVELOPMENT_PROGRAMS)

# Only list programs which need extra object files, or nonstandard linking
anchor_reads : genes.o utility.o
anchor_repeatness : files.o
assess_known_cn_events : assess_known_cn_events.gslo
bam_genotyper : utility.o
bridge_effect : genes.o utility.o
bridge_gene : genes.o utility.o
bridge_gene2 : genes.o utility.o
check_bridges : utility.o
check_bridges_result : utility.o
check_cand_cn : check_cand_cn.gslo
check_for_pseudogene : genes.o
check_mappability : files.o
compare_segmentation : compare_segmentation.gslo
compare_yoonha : utility.o
control_mview : control_mview.x11o
count_anchors : files.o
count_genes : genes.o utility.o
count_pseudogenes : files.o genes.o utility.o
debruijn : files.o genes.o utility.o
denovo_cn : denovo_cn.gslo
denovo_pseudogenes : files.o genes.o utility.o
encode : files.o
event_histogram : genes.o utility.o
family_count : utility.o
fastq_mapper : files.o utility.o
fastq_vt_mapper : files.o utility.o
find_bridge : utility.o
find_indel : utility.o
gene_view : genes.o utility.o
lookup_bridge : utility.o
mapper : files.o utility.o
mumdex_examples : utility.o
museq_primers : utility.o
optional_test : files.o
process_bisulfite_rna : files.o utility.o
primers : files.o
print_invariants : utility.o
rare : genes.o utility.o
repeatness : files.o
show_all_counts : files.o utility.o
show_counts : files.o utility.o
show_counts_special : files.o utility.o
show_mums : utility.o
show_pairs : utility.o
simulate_events : utility.o
simulate_events2 : utility.o
similar_bridge : utility.o
sssa : files.o
test_x11 : test_x11.x11o
transmission : genes.o utility.o
x11plot : x11plot.x11o



