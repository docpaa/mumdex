UTILITY_PROGRAMS := column_split \
	count_chars \
	count_words \
	echoe \
	even \
	even_columns \
	lines_at_pos \
	randomize_order \
	stats \
	table_columns \
	table_summary \
	test_beta \
	test_kd \
	test_lpe

PROGRAMS += $(UTILITY_PROGRAMS)
all : $(UTILITY_PROGRAMS)

# Only list programs which need extra object files, or nonstandard linking
test_beta: test_beta.gslo
