COMMON_LIBS=-L$(BIGHORNS)/lib -lfitslib $(NDIR)/slib/libbaselib.a $(NDIR)/slib/libmathlib.a  `root-config --libs`  -lcfitsio -lnova 
