# If you are using this Makefile standalone and fastjet-config is not
# in your path, edit this line to specify the full path
FASTJETCONFIG=fastjet-config
PREFIX=`$(FASTJETCONFIG) --prefix`
CXX=g++
CXXFLAGS= -O3 -Wall -g
install_script = ${SHELL} ../utils/install-sh
check_script = ../utils/check.sh

# global contrib-wide Makefile include may override some of the above
# variables (leading "-" means don't give an error if you can't find
# the file)
-include ../.Makefile.inc

#------------------------------------------------------------------------
# things that are specific to this contrib
NAME=JetFFMoments
SRCS=JetFFMoments.cc
INSTALLED_HEADERS=JetFFMoments.hh
#------------------------------------------------------------------------

# add fastjet libraries
CXXFLAGS+=$(shell ${FASTJETCONFIG} --cxxflags)
LDFLAGS+=-lm $(shell ${FASTJETCONFIG} --libs)

OBJS  = ${SRCS:.cc=.o}

install_HEADER  = $(install_script) -c -m 644
install_LIB     = $(install_script) -c -m 644
install_DIR     = $(install_script) -d

.PHONY: clean install examples

all: lib$(NAME).a

lib$(NAME).a: $(OBJS) 
	ar cru lib$(NAME).a $(OBJS)
	ranlib lib$(NAME).a

examples: example

example: all example.o
	$(CXX) -o example example.o -L. -l$(NAME) $(LDFLAGS)

check: examples
	@$(check_script) example ../data/Pythia-Zp2jets-lhc-pileup-1ev.dat || exit 1
	@echo "All tests successful"

examples: example

clean:
	rm -f *~ *.o

distclean: clean
	rm -f lib$(NAME).a example

install: all
	$(install_DIR) $(PREFIX)/include/fastjet/contrib
	$(install_HEADER) JetFFMoments.hh $(PREFIX)/include/fastjet/contrib
	$(install_DIR) $(PREFIX)/lib
	$(install_LIB) lib$(NAME).a $(PREFIX)/lib

depend:
	makedepend -Y --   -- $(SRCS) example.cc
# DO NOT DELETE

JetFFMoments.o: JetFFMoments.hh
example.o: JetFFMoments.hh
