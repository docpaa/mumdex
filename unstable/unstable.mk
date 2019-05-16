DEVELOPMENT_PROGRAMS := \
	add_abspos \
	afterburner \
	analyze_transmission \
	anchor_reads \
	anchor_mismatches \
	anchor_repeatness \
	annotate_candidates \
	annotate_repeats \
	assembler \
	assembly_index \
	assess_haha \
	assess_known_cn_events \
	averages \
	bam_genotyper \
	bin_events \
	bridge_counts \
	bridge_effect \
	bridge_gene2 \
	bridge_gene \
	bridge_properties \
	cancer_pop \
	candidate_stats \
	check_bridges \
	check_bridges_result \
	check_cand_cn \
	check_for_pseudogene \
	check_mappability \
	check_mumdex \
	colsplit \
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
	distinct_colors \
	dsDNAvt \
	encode \
	event_histogram \
	explore_coverage \
	family_count \
	family_figure \
	fastplot \
	fastq_kmer_counter \
	fastq_mapper \
	fastq_sequences \
	fastq_vt_mapper \
	filter_bridges \
	find_bridge \
	find_indel \
	find_microsatellite \
	find_repeats \
	gene_view \
	haha \
	hists \
	hist_compare \
	ihists \
	jackpot \
	karyotype \
	local_similarity \
	lookup_bridge \
	luria \
	make_reference \
	male_cn_stats \
	mapper \
	matrix \
	mismatches \
	missing_chromosome \
	model_transmission \
	mumdex_examples \
	mumdex_sequences \
	museq_primers \
	nlaIII_bsrsI \
	offset_beds \
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
	random_uniform_distances \
	rare \
	repeatness \
	repeats \
	sam_mapper \
	sd_errors \
	show_all_counts \
	show_counts \
	show_counts_special \
	show_mums \
	show_pairs \
	similar_bridge \
	simulate_cn \
	simulate_events2 \
	simulate_events \
	simulate_seating \
	skbr3_in_skn1 \
	smash_ip_fig \
	smooth_data \
	snp_report \
	subsample \
	table_stats \
	talk2018 \
	talk2018dec \
	test_hmm \
	test_mumdex \
	test_numerical \
	test_psplot \
	test_threads \
	test_x11 \
	transmission \
	transmission_old \
	unequal_bridges \
	validate_leukemia \
	validate_mumdex \
	venn \
	wg_smash_comparison \
	x11plot

PROGRAMS += $(DEVELOPMENT_PROGRAMS)
all : $(DEVELOPMENT_PROGRAMS)

# Only list programs which need extra object files, or nonstandard linking
anchor_reads : genes.o utility.o
afterburner: afterburner.gslo
anchor_repeatness : files.o
assess_known_cn_events : assess_known_cn_events.gslo
assess_haha : assess_haha.gslo
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
distinct_colors : distinct_colors.x11o
encode : files.o
event_histogram : genes.o utility.o
family_count : utility.o
fastq_mapper : files.o utility.o
fastq_vt_mapper : files.o utility.o
find_bridge : utility.o
find_indel : utility.o
gene_view : genes.o utility.o
haha: haha.gslo
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
sam_mapper : files.o utility.o
sd_errors : files.o
show_all_counts : files.o utility.o
show_counts : files.o utility.o
show_counts_special : files.o utility.o
show_mums : utility.o
show_pairs : utility.o
simulate_events : utility.o
simulate_events2 : utility.o
similar_bridge : utility.o
smash_ip_fig: files.o utility.o
sssa : files.o
test_x11 : test_x11.x11o
transmission : genes.o utility.o
transmission_old : genes.o utility.o
x11plot : x11plot.x11o



