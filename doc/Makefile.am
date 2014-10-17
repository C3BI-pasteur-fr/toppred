SUFFIXES = .pod .1

man_MANS = toppred.1

EXTRA_DIST = $(man_MANS) toppred.pod

PODARGS = --center='User Manuals' --release='Unix'

.pod.1:
	$(POD2MAN) $(PODARGS) $< > $@ && touch $@

## Maintainer parano checks
parano:
	(for p in *.pod; do podchecker --warnings --warnings $$p; done)
