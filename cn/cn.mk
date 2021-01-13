CN_PROGRAMS := \
	bad_cn_bins \
	combine_finebins \
	empirical_bins \
	extract_cn_segments \
	find_diploid_exons \
	finebin_cn \
	ggraph \
	heatmap \
	mumdex_cn \
	mumdex_finebin \
	pop_cn \
	recurrent_cn \
	smash \
	smash_qc \
	theoretical_bins

PROGRAMS += $(CN_PROGRAMS)
all : $(CN_PROGRAMS)

# Only list programs which need extra object files, or nonstandard linking
extract_cn_segments : extract_cn_segments.gslo genes.o
ggraph : ggraph.x11o genes.o
pop_cn: pop_cn.gslo genes.o
recurrent_cn: genes.o
smash : files.o utility.o
find_diploid_exons : genes.o

