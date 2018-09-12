CN_PROGRAMS := \
	bad_cn_bins \
	empirical_bins \
	extract_cn_segments \
	finebin_cn \
	ggraph \
	mumdex_cn \
	mumdex_finebin \
	pop_cn \
	smash \
	smash_qc

PROGRAMS += $(CN_PROGRAMS)
all : $(CN_PROGRAMS)

# Only list programs which need extra object files, or nonstandard linking
extract_cn_segments : extract_cn_segments.gslo genes.o
ggraph : ggraph.x11o genes.o
pop_cn: pop_cn.gslo genes.o
smash : files.o utility.o

