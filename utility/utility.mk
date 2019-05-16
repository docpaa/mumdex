UTILITY_PROGRAMS := column_split \
	even \
	even_columns \
	lines_at_pos \
	randomize_order \
	stats

PROGRAMS += $(UTILITY_PROGRAMS)
all : $(UTILITY_PROGRAMS)

# Only list programs which need extra object files, or nonstandard linking
