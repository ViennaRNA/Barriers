manpages_verbose = $(manpages_verbose_@AM_V@)
manpages_verbose_ = $(manpages_verbose_@AM_DEFAULT_V@)
manpages_verbose_0 = @echo "  MAN      $@";

dist_man_MANS = barriers.1

SUFFIXES = .1 .ggo

vpath %.ggo = ../src

if BARRIERS_AM_SWITCH_HAS_HELP2MAN

.ggo.1 :
	$(manpages_verbose)./help2manWrapper.sh \
                "$(HELP2MAN)" $< | perl -pe 's/\\e\\e/\\e/g' > $@

endif

EXTRA_DIST = \
    include \
    barriers.texinfo

# remove all generated man pages
maintainer-clean-local:
	-rm -f $(dist_man_MANS)

