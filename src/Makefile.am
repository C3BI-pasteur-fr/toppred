
AM_CPPFLAGS = -DDATADIR=\"$(pkgdatadir)\" $(OS_CPPFLAGS)

bin_PROGRAMS = toppred

if USE_GD_LIB_SRC
extra_src=graph.c
extra_header=graph.h
else
extra_src=
extra_header=
endif

EXTRA_DIST = graph.c graph.h

toppred_SOURCES = error.c main.c usage.c profile.c seq-reader.c	loop.c \
	charge.c topology.c topoprint.c	output.c mloutput.c $(extra_src)

noinst_HEADERS = error.h main.h usage.h profile.h seq-reader.h loop.h \
	params.h charge.h topology.h topoprint.h output.h mloutput.h \
	$(extra_header)

## Maintainer parano check
LINTDEFS = $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS)

parano:
	$(LINT) $(LINTFLAGS) $(LINTDEFS) $(toppred_SOURCES)
