# If you are using this Makefile standalone and fastjet-config is not
# in your path, edit this line to specify the full path
FASTJETCONFIG=fastjet-config
PREFIX=`$(FASTJETCONFIG) --prefix`
CXX=g++
CXXFLAGS= -O3 -Wall -g
install_script = $(SHELL) ../utils/install-sh
check_script = ../utils/check.sh

# global contrib-wide Makefile include may override some of the above
# variables (leading "-" means don't give an error if you can't find
# the file)
-include ../.Makefile.inc

#------------------------------------------------------------------------
# things that are specific to this contrib
NAME=RecursiveTools
SRCS=Recluster.cc RecursiveSymmetryCutBase.cc ModifiedMassDropTagger.cc SoftDrop.cc IteratedSoftDrop.cc RecursiveSoftDrop.cc BottomUpSoftDrop.cc
EXAMPLES=example_mmdt example_mmdt_sub example_mmdt_ee example_recluster example_softdrop example_recursive_softdrop example_advanced_usage example_isd example_bottomup_softdrop
INSTALLED_HEADERS=Recluster.hh RecursiveSymmetryCutBase.hh ModifiedMassDropTagger.hh SoftDrop.hh IteratedSoftDrop.hh RecursiveSoftDrop.hh BottomUpSoftDrop.hh
#------------------------------------------------------------------------

CXXFLAGS+= $(shell $(FASTJETCONFIG) --cxxflags)
LDFLAGS += -lm $(shell $(FASTJETCONFIG) --libs)

OBJS  = $(SRCS:.cc=.o)
EXAMPLES_SRCS  = $(EXAMPLES:=.cc)

install_HEADER  = $(install_script) -c -m 644
install_LIB     = $(install_script) -c -m 644
install_DIR     = $(install_script) -d
install_DATA    = $(install_script) -c -m 644
install_PROGRAM = $(install_script) -c -s
install_SCRIPT  = $(install_script) -c

.PHONY: clean distclean examples check install

# compilation of the code (default target)
all: lib$(NAME).a

lib$(NAME).a: $(OBJS) 
	ar cru lib$(NAME).a $(OBJS)
	ranlib lib$(NAME).a

# building the examples
examples: $(EXAMPLES)

# the following construct makes it possible to automatically build
# each of the examples listed in $EXAMPLES
$(EXAMPLES): % : %.o all
	$(CXX) -o $@ $< -L. -l$(NAME) $(LDFLAGS)

# check that everything went fine
check: examples
	@$(check_script) example_mmdt ../data/single-event.dat || exit 1
	@$(check_script) example_mmdt_sub  ../data/Pythia-Zp2jets-lhc-pileup-1ev.dat || exit 1
	@$(check_script) example_mmdt_ee  ../data/single-ee-event.dat || exit 1
	@$(check_script) example_recluster ../data/single-event.dat || exit 1
	@$(check_script) example_softdrop  ../data/single-event.dat || exit 1
	@$(check_script) example_advanced_usage  ../data/single-event.dat || exit 1
	@$(check_script) example_recursive_softdrop  ../data/single-event.dat || exit 1
	@$(check_script) example_bottomup_softdrop  ../data/single-event.dat || exit 1
	@$(check_script) example_isd  ../data/single-event.dat || exit 1
#	@for prog in $(EXAMPLES); do\
#	  $(check_script) $${prog} ../data/single-event.dat || exit 1; \
#	done
	@echo "All tests successful"

# cleaning the directory
clean:
	rm -f *~ *.o

distclean: clean
	rm -f lib$(NAME).a $(EXAMPLES)

# install things in PREFIX/...
install: all
	$(install_DIR) $(PREFIX)/include/fastjet/contrib
	for header in $(INSTALLED_HEADERS); do\
	  $(install_HEADER) $$header $(PREFIX)/include/fastjet/contrib/;\
	done
	$(install_DIR) $(PREFIX)/lib
	$(install_LIB) lib$(NAME).a $(PREFIX)/lib

depend:
	makedepend -Y --   -- $(SRCS) $(EXAMPLES_SRCS)
# DO NOT DELETE

Recluster.o: Recluster.hh
RecursiveSymmetryCutBase.o: RecursiveSymmetryCutBase.hh Recluster.hh
ModifiedMassDropTagger.o: ModifiedMassDropTagger.hh
ModifiedMassDropTagger.o: RecursiveSymmetryCutBase.hh Recluster.hh
SoftDrop.o: SoftDrop.hh RecursiveSymmetryCutBase.hh Recluster.hh
IteratedSoftDrop.o: IteratedSoftDrop.hh
RecursiveSoftDrop.o: RecursiveSoftDrop.hh Recluster.hh SoftDrop.hh
RecursiveSoftDrop.o: RecursiveSymmetryCutBase.hh
BottomUpSoftDrop.o: BottomUpSoftDrop.hh
example_mmdt.o: ModifiedMassDropTagger.hh RecursiveSymmetryCutBase.hh
example_mmdt.o: Recluster.hh
example_mmdt_sub.o: ModifiedMassDropTagger.hh RecursiveSymmetryCutBase.hh
example_mmdt_sub.o: Recluster.hh
example_recluster.o: Recluster.hh
example_softdrop.o: SoftDrop.hh RecursiveSymmetryCutBase.hh Recluster.hh
example_recursive_softdrop.o: RecursiveSoftDrop.hh Recluster.hh SoftDrop.hh
example_recursive_softdrop.o: RecursiveSymmetryCutBase.hh
example_bottomup_softdrop.o: BottomUpSoftDrop.hh
example_advanced_usage.o: SoftDrop.hh RecursiveSymmetryCutBase.hh
example_advanced_usage.o: Recluster.hh
example_isd.o: IteratedSoftDrop.hh
