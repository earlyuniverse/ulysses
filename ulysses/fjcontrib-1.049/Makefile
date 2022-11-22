# -*- Makefile -*-
# This Makefile was generated automatically by the 'configure' script using
#    ./configure --fastjet-config=/Users/jessicaturner/Documents/GitHub/ulysses/ulysses/local/bin/fastjet-config CXXFLAGS=-fPIC
# If you edit the processed version, any changes will be lost after the next configure
# (for developers, make sure you only edit Makefile.in)

# installation setup
SUBDIRS=JetFFMoments ClusteringVetoPlugin LundPlane ConstituentSubtractor Centauro SoftKiller EnergyCorrelator JetsWithoutJets ScJet Nsubjettiness QCDAwarePlugin GenericSubtractor ValenciaPlugin FlavorCone VariableR RecursiveTools SubjetCounting JetCleanser
SUBDIRS.all=$(SUBDIRS:=.all)

# these will be overriden if the user specifies CXX or CXXFLAGS with configure
# they will also be overriden by definitions in subsiduary Makefiles
CXX=g++
CXXFLAGS=-O2 -Wall -g

# get any variables defined in the contrib-wide include
-include .Makefile.inc

# make invocations in subdirectories have an environment variable
# set so that they can alter their behaviour if need be 
# (e.g. for make check)
SUBMAKE= $(MAKE) FJCONTRIB_SUBMAKE=1

.PHONY: $(SUBDIRS) $(SUBDIRS.all) clean distclean check check_init install examples

all: $(SUBDIRS.all)

install: all $(SUBDIRS)

examples clean: $(SUBDIRS)

distclean: $(SUBDIRS) fragile-shared-distclean

check: check_init $(SUBDIRS)
	@echo ""
	@cat test_summary.tmp
	@printf "\n%d out of %d tests passed\n\n" `grep "Success" test_summary.tmp | wc -l` `grep "^  " test_summary.tmp | wc -l`
	@rm -f test_summary.tmp

check_init:
	@echo "Summary of tests" >  test_summary.tmp
	@echo "----------------" >> test_summary.tmp

# distclean removes the Makefile, but leaves in config.log
distclean:
	rm -f Makefile
	rm -f .Makefile.inc

# dirty hack to provide a shared library to CMS; this is extremely fragile
# and will be hopefully replaced with a more robust solution at some
# point in the future; note, that you must still do make install
# in order for all the headers to be installed
fragile-shared: libfastjetcontribfragile.dylib

fragile_SHARED_SRC_LIST=JetFFMoments/JetFFMoments.cc ClusteringVetoPlugin/ClusteringVetoPlugin.cc LundPlane/LundGenerator.cc LundPlane/LundWithSecondary.cc LundPlane/SecondaryLund.cc LundPlane/RecursiveLundEEGenerator.cc ConstituentSubtractor/ConstituentSubtractor.cc ConstituentSubtractor/IterativeConstituentSubtractor.cc ConstituentSubtractor/RescalingClasses.cc Centauro/Centauro.cc SoftKiller/SoftKiller.cc EnergyCorrelator/EnergyCorrelator.cc JetsWithoutJets/JetsWithoutJets.cc JetsWithoutJets/EventStorage.cc ScJet/ScJet.cc Nsubjettiness/Nsubjettiness.cc Nsubjettiness/Njettiness.cc Nsubjettiness/NjettinessPlugin.cc Nsubjettiness/XConePlugin.cc Nsubjettiness/MeasureDefinition.cc Nsubjettiness/ExtraRecombiners.cc Nsubjettiness/AxesDefinition.cc Nsubjettiness/TauComponents.cc QCDAwarePlugin/QCDAwarePlugin.cc GenericSubtractor/GenericSubtractor.cc ValenciaPlugin/ValenciaPlugin.cc FlavorCone/FlavorCone.cc VariableR/VariableRPlugin.cc RecursiveTools/Recluster.cc RecursiveTools/RecursiveSymmetryCutBase.cc RecursiveTools/ModifiedMassDropTagger.cc RecursiveTools/SoftDrop.cc RecursiveTools/IteratedSoftDrop.cc RecursiveTools/RecursiveSoftDrop.cc RecursiveTools/BottomUpSoftDrop.cc SubjetCounting/SubjetCounting.cc JetCleanser/JetCleanser.cc
libfastjetcontribfragile.dylib: $(fragile_SHARED_SRC_LIST)
	$(CXX) -dynamiclib -fPIC -DPIC $(CXXFLAGS) `$(FASTJETCONFIG) --cxxflags --libs` $(fragile_SHARED_SRC_LIST) -o libfastjetcontribfragile.dylib

fragile-shared-install: fragile-shared
	utils/install-sh -c -m 755 libfastjetcontribfragile.dylib $(PREFIX)/lib
	install_name_tool -id $(PREFIX)/lib/libfastjetcontribfragile.dylib  $(PREFIX)/lib/libfastjetcontribfragile.dylib 

fragile-shared-distclean:
	rm -f libfastjetcontribfragile.dylib

$(SUBDIRS):
	+$(SUBMAKE) -C $@ $(MAKECMDGOALS)

$(SUBDIRS.all):
	+$(SUBMAKE) -C $(basename $@)
