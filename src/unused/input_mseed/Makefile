MODULE_NAME = input_mseed

MODDIR    = $(SRCDIR)/cs/modules/$(MODULE_NAME)
MODULE    = $(MODDIR)/mod_$(MODULE_NAME).cc
OBJS      = $(OBJDIR)/mod_$(MODULE_NAME).o
LIB_v1.0  = libmod_$(MODULE_NAME).so.1.0

INCS =  -I"$(SRCDIR)/cs/geolib"  -I"$(SRCDIR)/cs/system" -I"$(MODDIR)"

ALL_FLAGS = $(INCS) $(COMMON_FLAGS) -fPIC


OBJ_MSEED  = $(OBJDIR)/fileutils.o $(OBJDIR)/genutils.o $(OBJDIR)/gswap.o $(OBJDIR)/lookup.o $(OBJDIR)/msrutils.o $(OBJDIR)/pack.o $(OBJDIR)/unpack.o $(OBJDIR)/packdata.o $(OBJDIR)/unpackdata.o $(OBJDIR)/traceutils.o $(OBJDIR)/lmplatform.o $(OBJDIR)/MSeedReader.o $(OBJDIR)/FlexHeader.o

default: $(OBJ_MSEED) $(OBJS) $(LIBDIR)/$(LIB_v1.0)

clean:
	rm -f $(OBJ_MSEED) $(OBJS)

bleach: clean
	rm -f $(LIBDIR)/$(LIB_v1.0)

$(OBJDIR)/mod_$(MODULE_NAME).o: $(MODULE)
	$(CPP) -c $(ALL_FLAGS) $(MODULE) -o $(OBJDIR)/mod_$(MODULE_NAME).o

$(LIBDIR)/$(LIB_v1.0): $(OBJ_MSEED) $(OBJS)
	$(CPP) -shared -Wl,-$(SONAME),$(LIB_v1.0) -o $(LIBDIR)/$(LIB_v1.0) $(OBJ_MSEED) $(OBJS) -L$(LIBDIR) -lc -lgeolib -lcseis_system


$(OBJDIR)/fileutils.o:	$(MODDIR)/fileutils.c
			$(CC) $(ALL_FLAGS) $(MODDIR)/fileutils.c -c -o $(OBJDIR)/fileutils.o

$(OBJDIR)/genutils.o:	$(MODDIR)/genutils.c
			$(CC) $(ALL_FLAGS) $(MODDIR)/genutils.c -c -o $(OBJDIR)/genutils.o

$(OBJDIR)/lmplatform.o:	$(MODDIR)/lmplatform.c
			$(CC) $(ALL_FLAGS) $(MODDIR)/lmplatform.c -c -o $(OBJDIR)/lmplatform.o

$(OBJDIR)/lookup.o:	$(MODDIR)/lookup.c
			$(CC) $(ALL_FLAGS) $(MODDIR)/lookup.c -c -o $(OBJDIR)/lookup.o

$(OBJDIR)/msrutils.o:	$(MODDIR)/msrutils.c
			$(CC) $(ALL_FLAGS) $(MODDIR)/msrutils.c -c -o $(OBJDIR)/msrutils.o

$(OBJDIR)/pack.o:	$(MODDIR)/pack.c
			$(CC) $(ALL_FLAGS) $(MODDIR)/pack.c -c -o $(OBJDIR)/pack.o

$(OBJDIR)/gswap.o:	$(MODDIR)/gswap.c
			$(CC) $(ALL_FLAGS) $(MODDIR)/gswap.c -c -o $(OBJDIR)/gswap.o

$(OBJDIR)/packdata.o:	$(MODDIR)/packdata.c
			$(CC) $(ALL_FLAGS) $(MODDIR)/packdata.c -c -o $(OBJDIR)/packdata.o

$(OBJDIR)/traceutils.o:	$(MODDIR)/traceutils.c
			$(CC) $(ALL_FLAGS) $(MODDIR)/traceutils.c -c -o $(OBJDIR)/traceutils.o

$(OBJDIR)/unpack.o:	$(MODDIR)/unpack.c
			$(CC) $(ALL_FLAGS) $(MODDIR)/unpack.c -c -o $(OBJDIR)/unpack.o

$(OBJDIR)/unpackdata.o:	$(MODDIR)/unpackdata.c
			$(CC) $(ALL_FLAGS) $(MODDIR)/unpackdata.c -c -o $(OBJDIR)/unpackdata.o

$(OBJDIR)/MSeedReader.o:	$(MODDIR)/MSeedReader.cc $(MODDIR)/MSeedReader.h
			$(CPP) $(ALL_FLAGS) $(MODDIR)/MSeedReader.cc -c -o $(OBJDIR)/MSeedReader.o

$(OBJDIR)/FlexHeader.o:	$(MODDIR)/FlexHeader.cc $(MODDIR)/FlexHeader.h
			$(CPP) $(ALL_FLAGS) $(MODDIR)/FlexHeader.cc -c -o $(OBJDIR)/FlexHeader.o
