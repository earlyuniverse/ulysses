# If you are using this Makefile standalone and fastjet-config is not
# in your path, edit this line to specify the full path
FASTJETCONFIG=fastjet-config
PREFIX=`$(FASTJETCONFIG) --prefix`
CXX=g++
CXXFLAGS= -O3 -Wall -g -fPIC
install_script = $(SHELL) ../utils/install-sh
check_script = ../utils/check.sh

# global contrib-wide Makefile include may override some of the above
# variables (leading "-" means don't give an error if you can't find
# the file)
-include ../.Makefile.inc

#------------------------------------------------------------------------
# things that are specific to this contrib
NAME=ConstituentSubtractor
SRCS=ConstituentSubtractor.cc IterativeConstituentSubtractor.cc RescalingClasses.cc
EXAMPLES=example_jet_by_jet example_event_wide example_iterative example_background_rescaling example_whole_event_using_charged_info
INSTALLED_HEADERS=ConstituentSubtractor.hh IterativeConstituentSubtractor.hh RescalingClasses.hh
#------------------------------------------------------------------------

CXXFLAGS+= $(shell $(FASTJETCONFIG) --cxxflags)
LDFLAGS += -lm $(shell $(FASTJETCONFIG) --libs)
#CXXFLAGS+= $(shell root-config --cflags)
#LDFLAGS+= -lm $(shell root-config --libs)

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
	@for prog in $(EXAMPLES); do\
	  $(check_script) $${prog} ../data/Pythia-Zp2jets-lhc-pileup-1ev.dat || echo "failed"; \
	done
#	@echo "All tests successful"
#	  $(check_script) $${prog} ../data/Pythia-Zp2jets-lhc-pileup-1ev.dat || exit 1; \

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
