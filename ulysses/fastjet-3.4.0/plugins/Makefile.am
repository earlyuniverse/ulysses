SUBDIRS = 

if BUILD_PLUGIN_SISCONE
SUBDIRS += SISCone
endif
if BUILD_PLUGIN_CDFCONES
SUBDIRS += CDFCones
endif
if BUILD_PLUGIN_PXCONE
SUBDIRS += PxCone
endif
if BUILD_PLUGIN_D0RUNIICONE
SUBDIRS += D0RunIICone
endif
if BUILD_PLUGIN_NESTEDDEFS
SUBDIRS += NestedDefs
endif
if BUILD_PLUGIN_TRACKJET
SUBDIRS += TrackJet
endif
if BUILD_PLUGIN_ATLASCONE
SUBDIRS += ATLASCone
endif
if BUILD_PLUGIN_EECAMBRIDGE
SUBDIRS += EECambridge
endif
if BUILD_PLUGIN_JADE
SUBDIRS += Jade
endif
if BUILD_PLUGIN_CMSITERATIVECONE
SUBDIRS += CMSIterativeCone
endif
if BUILD_PLUGIN_D0RUNICONE
SUBDIRS += D0RunICone
endif
if BUILD_PLUGIN_GRIDJET
SUBDIRS += GridJet
endif
# end plugin section -- do not modify this line

SUBDIRS += .

if MONOLITHIC_PLUGINS
lib_LTLIBRARIES = libfastjetplugins.la
libfastjetplugins_la_SOURCES=
nodist_EXTRA_libfastjetplugins_la_SOURCES = dummy.cc
libfastjetplugins_la_LIBADD = ../src/libfastjet.la
if BUILD_PLUGIN_SISCONE
libfastjetplugins_la_LIBADD += SISCone/libSISConePlugin.la
endif
if BUILD_PLUGIN_CDFCONES
libfastjetplugins_la_LIBADD += CDFCones/libCDFConesPlugin.la
endif
if BUILD_PLUGIN_PXCONE
libfastjetplugins_la_LIBADD += PxCone/libPxConePlugin.la
endif
if BUILD_PLUGIN_D0RUNIICONE
libfastjetplugins_la_LIBADD += D0RunIICone/libD0RunIIConePlugin.la
endif
if BUILD_PLUGIN_NESTEDDEFS
libfastjetplugins_la_LIBADD += NestedDefs/libNestedDefsPlugin.la
endif
if BUILD_PLUGIN_TRACKJET
libfastjetplugins_la_LIBADD += TrackJet/libTrackJetPlugin.la
endif
if BUILD_PLUGIN_ATLASCONE
libfastjetplugins_la_LIBADD += ATLASCone/libATLASConePlugin.la
endif
if BUILD_PLUGIN_EECAMBRIDGE
libfastjetplugins_la_LIBADD += EECambridge/libEECambridgePlugin.la
endif
if BUILD_PLUGIN_JADE
libfastjetplugins_la_LIBADD += Jade/libJadePlugin.la
endif
if BUILD_PLUGIN_CMSITERATIVECONE
libfastjetplugins_la_LIBADD += CMSIterativeCone/libCMSIterativeConePlugin.la
endif
if BUILD_PLUGIN_D0RUNICONE
libfastjetplugins_la_LIBADD += D0RunICone/libD0RunIConePlugin.la
endif
if BUILD_PLUGIN_GRIDJET
libfastjetplugins_la_LIBADD += GridJet/libGridJetPlugin.la
endif
# end monolithic plugin section -- do not modify this line
endif

DIST_SUBDIRS = GridJet D0RunICone CMSIterativeCone Jade EECambridge ATLASCone TrackJet SISCone CDFCones PxCone D0RunIICone NestedDefs

EXTRA_DIST = makefile.static
