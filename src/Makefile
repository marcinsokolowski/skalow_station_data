# include flags.mak
# include $(BIGHORNS)/software/analysis/make.cfg/flags.mak

# pre-processor :
COMMON_LIBS=-L$(BIGHORNS)/lib -lfitslib -lmsfitslib  -lcfitsio -lnova 
INCLUDES=-D_UNIX -I$(SRCDIR)/cmn/mathlib/ -I$(SRCDIR)/cmn/baselib/ -I$(BIGHORNS)/software/analysis/fitslib/

all :
	g++ read_binary_station_beam_test_order1_2pol.c $(OPT) -o read_binary_station_beam_test_order1_2pol -lm
	g++ read_binary_station_beam_test_order1_2pol_totpower.c  $(OPT) -o read_binary_station_beam_test_order1_2pol_totpower -lm
	g++ test_pattern_generated_data.c $(OPT) -o test_pattern_generated_data -lm
	g++ skalow_spectrometer.cpp -D_UNIX $(OPT) -o skalow_spectrometer -L$(BIGHORNS)/lib $(COMMON_LIBS) -lm -lfftw3 -lmsfitslib
#	g++ skalow_stitch_data.cpp -D_UNIX $(OPT) -o skalow_stitch_data -L$(BIGHORNS)/lib -lfitslib $(COMMON_LIBS) -lm
#	g++ skalow_stitch_data_transposed.cpp -D_UNIX $(OPT) -o skalow_stitch_data_transposed -L$(BIGHORNS)/lib -lfitslib $(COMMON_LIBS) -lm
	
	cp read_binary_station_beam_test_order1_2pol $(NDIR)/bin/
	cp read_binary_station_beam_test_order1_2pol_totpower $(NDIR)/bin/
	cp test_pattern_generated_data $(NDIR)/bin/
	cp skalow_spectrometer $(NDIR)/bin/
#	cp skalow_stitch_data $(NDIR)/bin/
#	cp skalow_stitch_data_transposed $(NDIR)/bin/
	

clean :
	rm -f read_binary_station_beam_test_order1_2pol test_pattern_generated_data read_binary_station_beam_test_order1_2pol_totpower *.o skalow_spectrometer skalow_stitch_data skalow_stitch_data_transposed
	