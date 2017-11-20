UTILITY_PROGRAMS := even_columns randomize_order
PROGRAMS += $(UTILITY_PROGRAMS)
all : $(UTILITY_PROGRAMS)

# Only list programs which need extra object files, or nonstandard linking
