# for the time being allow these things to be distributed but without
# any particular help to go with them? This should change at some point

SUBDIRS = data python root
EXTRA_DIST = Makefile.alt ktjet_example.cc ktjet_timing.cc CmdLine.hh\
	makefile.static
#	test-example.sh test-areas.sh makefile.static

noinst_PROGRAMS = fastjet_example fastjet_areas fastjet_timing_plugins\
 01-basic  \
 02-jetdef \
 04-constituents \
 05-eplus_eminus \
 06-area \
 07-subtraction \
 08-selector \
 09-user_info \
 10-subjets \
 11-filter \
 12-boosted_higgs \
 13-boosted_top \
 14-groomers

check_PROGRAMS = fastjet_example fastjet_areas

# basic usage example
01_basic_SOURCES = 01-basic.cc
01_basic_SHORTNAME = 01
01_basic_CXXFLAGS = $(AM_CXXFLAGS) $(CGAL_CPPFLAGS) -I$(srcdir)/../include
01_basic_LDADD = ../src/libfastjet.la $(CGAL_LIBS)

# further illustrates the usage of jet definitions
02_jetdef_SOURCES = 02-jetdef.cc
02_jetdef_SHORTNAME = 02
02_jetdef_CXXFLAGS = $(AM_CXXFLAGS) $(CGAL_CPPFLAGS) -I$(srcdir)/../include
02_jetdef_LDADD = ../src/libfastjet.la $(CGAL_LIBS)

# illustrates the usage of plugins
if BUILD_PLUGIN_SISCONE
noinst_PROGRAMS += 03-plugin
03_plugin_SOURCES = 03-plugin.cc
#03_plugin_SHORTNAME = 03
03_plugin_CXXFLAGS = $(AM_CXXFLAGS) $(CGAL_CPPFLAGS) -I$(srcdir)/../include
03_plugin_LDADD = ../src/libfastjet.la $(CGAL_LIBS)
03_plugin_CXXFLAGS += -I$(srcdir)/../plugins/SISCone -I$(srcdir)/../plugins/SISCone/siscone
03_plugin_LDADD += ../plugins/SISCone/libSISConePlugin.la ../plugins/SISCone/siscone/siscone/libsiscone.la ../plugins/SISCone/siscone/siscone/spherical/libsiscone_spherical.la
else
EXTRA_DIST += 03-plugin.cc
endif

# illustrate the access to the constituents
04_constituents_SOURCES = 04-constituents.cc
04_constituents_SHORTNAME = 04
04_constituents_CXXFLAGS = $(AM_CXXFLAGS) $(CGAL_CPPFLAGS) -I$(srcdir)/../include
04_constituents_LDADD = ../src/libfastjet.la $(CGAL_LIBS)

# illustrate clustering in e+e- collisions
05_eplus_eminus_SOURCES = 05-eplus_eminus.cc
05_eplus_eminus_SHORTNAME = 05
05_eplus_eminus_CXXFLAGS = $(AM_CXXFLAGS) $(CGAL_CPPFLAGS) -I$(srcdir)/../include
05_eplus_eminus_LDADD = ../src/libfastjet.la $(CGAL_LIBS)

# illustrate the computation of jet areas
06_area_SOURCES = 06-area.cc
06_area_SHORTNAME = 06
06_area_CXXFLAGS = $(AM_CXXFLAGS) $(CGAL_CPPFLAGS) -I$(srcdir)/../include
06_area_LDADD = ../src/libfastjet.la $(CGAL_LIBS)

# illustrate the background subtraction using jet areas
07_subtraction_SOURCES = 07-subtraction.cc
07_subtraction_SHORTNAME = 07
07_subtraction_CXXFLAGS = $(AM_CXXFLAGS) $(CGAL_CPPFLAGS) -I$(srcdir)/../include -I$(srcdir)/../tools
07_subtraction_LDADD = ../tools/libfastjettools.la ../src/libfastjet.la $(CGAL_LIBS)

# illustrate the use of Selector's
08_selector_SOURCES = 08-selector.cc
08_selector_SHORTNAME = 08
08_selector_CXXFLAGS = $(AM_CXXFLAGS) $(CGAL_CPPFLAGS) -I$(srcdir)/../include 
08_selector_LDADD =  ../src/libfastjet.la $(CGAL_LIBS)

# illustrate how to associate user information to the jets
09_user_info_SOURCES = 09-user_info.cc
09_user_info_SHORTNAME = 09
09_user_info_CXXFLAGS = $(AM_CXXFLAGS) $(CGAL_CPPFLAGS) -I$(srcdir)/../include -I$(srcdir)/../tools
09_user_info_LDADD = ../tools/libfastjettools.la ../src/libfastjet.la $(CGAL_LIBS)

# illustrate how to access subjets
10_subjets_SOURCES = 10-subjets.cc
10_subjets_SHORTNAME = 10
10_subjets_CXXFLAGS = $(AM_CXXFLAGS) $(CGAL_CPPFLAGS) -I$(srcdir)/../include -I$(srcdir)/../tools
10_subjets_LDADD = ../src/libfastjet.la $(CGAL_LIBS)

# illustrate the use of Filter
11_filter_SOURCES = 11-filter.cc
11_filter_SHORTNAME = 11
11_filter_CXXFLAGS = $(AM_CXXFLAGS) $(CGAL_CPPFLAGS) -I$(srcdir)/../include -I$(srcdir)/../tools
11_filter_LDADD = ../tools/libfastjettools.la ../src/libfastjet.la  $(CGAL_LIBS)

# illustrate how to tag boosted higgs
12_boosted_higgs_SOURCES = 12-boosted_higgs.cc
12_boosted_higgs_SHORTNAME = 12
12_boosted_higgs_CXXFLAGS = $(AM_CXXFLAGS) $(CGAL_CPPFLAGS) -I$(srcdir)/../include -I$(srcdir)/../tools
12_boosted_higgs_LDADD = ../tools/libfastjettools.la ../src/libfastjet.la $(CGAL_LIBS)

# illustrate how to tag boosted top
13_boosted_top_SOURCES = 13-boosted_top.cc
13_boosted_top_SHORTNAME = 13
13_boosted_top_CXXFLAGS = $(AM_CXXFLAGS) $(CGAL_CPPFLAGS) -I$(srcdir)/../include -I$(srcdir)/../tools
13_boosted_top_LDADD = ../tools/libfastjettools.la ../src/libfastjet.la  $(CGAL_LIBS)

# illustrate how to use various groomers in a unified manner
14_groomers_SOURCES = 14-groomers.cc
14_groomers_SHORTNAME = 14
14_groomers_CXXFLAGS = $(AM_CXXFLAGS) $(CGAL_CPPFLAGS) -I$(srcdir)/../include -I$(srcdir)/../tools
14_groomers_LDADD = ../tools/libfastjettools.la ../src/libfastjet.la  $(CGAL_LIBS)


# backward-compatible + make-check-related things
fastjet_example_SOURCES = fastjet_example.cc
fastjet_example_CXXFLAGS = $(AM_CXXFLAGS) $(CGAL_CPPFLAGS) -I$(srcdir)/../include
fastjet_example_LDADD = ../src/libfastjet.la $(CGAL_LIBS)

fastjet_areas_SOURCES = fastjet_areas.cc
fastjet_areas_CXXFLAGS = $(AM_CXXFLAGS) $(CGAL_CPPFLAGS) -I$(srcdir)/../include
fastjet_areas_LDADD = ../src/libfastjet.la $(CGAL_LIBS)
if BUILD_PLUGIN_SISCONE
fastjet_areas_CXXFLAGS += -I$(srcdir)/../plugins/SISCone -I$(srcdir)/../plugins/SISCone/siscone
fastjet_areas_LDADD += ../plugins/SISCone/libSISConePlugin.la ../plugins/SISCone/siscone/siscone/libsiscone.la ../plugins/SISCone/siscone/siscone/spherical/libsiscone_spherical.la
endif

fastjet_timing_plugins_SOURCES = fastjet_timing_plugins.cc CmdLine.cc
fastjet_timing_plugins_CXXFLAGS = $(AM_CXXFLAGS) $(CGAL_CPPFLAGS) -I$(srcdir)/../include
fastjet_timing_plugins_LDADD = ../src/libfastjet.la $(CGAL_LIBS)
fastjet_timing_plugins_CXXFLAGS += -I$(srcdir)/../tools
fastjet_timing_plugins_LDADD += ../tools/libfastjettools.la 
if BUILD_PLUGIN_SISCONE
fastjet_timing_plugins_CXXFLAGS += -I$(srcdir)/../plugins/SISCone -I$(srcdir)/../plugins/SISCone/siscone
fastjet_timing_plugins_LDADD += ../plugins/SISCone/libSISConePlugin.la ../plugins/SISCone/siscone/siscone/libsiscone.la ../plugins/SISCone/siscone/siscone/spherical/libsiscone_spherical.la
endif
if BUILD_PLUGIN_CDFCONES
fastjet_timing_plugins_CXXFLAGS += -I$(srcdir)/../plugins/CDFCones
fastjet_timing_plugins_LDADD += ../plugins/CDFCones/libCDFConesPlugin.la ../plugins/CDFCones/CDFcode/libCDFcode.la
endif
if BUILD_PLUGIN_PXCONE
fastjet_timing_plugins_CXXFLAGS += -I$(srcdir)/../plugins/PxCone
fastjet_timing_plugins_LDADD += ../plugins/PxCone/libPxConePlugin.la $(FLIBS)
endif
if BUILD_PLUGIN_D0RUNIICONE
fastjet_timing_plugins_CXXFLAGS += -I$(srcdir)/../plugins/D0RunIICone
fastjet_timing_plugins_LDADD += ../plugins/D0RunIICone/libD0RunIIConePlugin.la $(FLIBS)
endif
if BUILD_PLUGIN_TRACKJET
fastjet_timing_plugins_CXXFLAGS += -I$(srcdir)/../plugins/TrackJet
fastjet_timing_plugins_LDADD += ../plugins/TrackJet/libTrackJetPlugin.la
endif
if BUILD_PLUGIN_ATLASCONE
fastjet_timing_plugins_CXXFLAGS += -I$(srcdir)/../plugins/ATLASCone
fastjet_timing_plugins_LDADD += ../plugins/ATLASCone/libATLASConePlugin.la
endif
if BUILD_PLUGIN_EECAMBRIDGE
fastjet_timing_plugins_CXXFLAGS += -I$(srcdir)/../plugins/EECambridge
fastjet_timing_plugins_LDADD += ../plugins/EECambridge/libEECambridgePlugin.la
endif
if BUILD_PLUGIN_JADE
fastjet_timing_plugins_CXXFLAGS += -I$(srcdir)/../plugins/Jade
fastjet_timing_plugins_LDADD += ../plugins/Jade/libJadePlugin.la
endif
if BUILD_PLUGIN_CMSITERATIVECONE
fastjet_timing_plugins_CXXFLAGS += -I$(srcdir)/../plugins/CMSIterativeCone
fastjet_timing_plugins_LDADD += ../plugins/CMSIterativeCone/libCMSIterativeConePlugin.la
endif
if BUILD_PLUGIN_D0RUNICONE
fastjet_timing_plugins_CXXFLAGS += -I${srcdir}/../plugins/D0RunICone
fastjet_timing_plugins_LDADD += ../plugins/D0RunICone/libD0RunIConePlugin.la
endif
if BUILD_PLUGIN_GRIDJET
fastjet_timing_plugins_CXXFLAGS += -I${srcdir}/../plugins/GridJet
fastjet_timing_plugins_LDADD += ../plugins/GridJet/libGridJetPlugin.la
endif

# add entries for other plugins here

#TESTS_ENVIRONMENT = $(SHELL)
#TESTS = test-example.sh test-areas.sh
