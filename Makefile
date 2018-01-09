CC=gcc
CXX=g++
RM=rm -f
CPPFLAGS=-03
INCLUDES=-Iinc
LDFLAGS=
LDLIBS=-lglpk

SRCDIR=src
OBJDIR=obj
INCDIR=inc
BINDIR=bin
TESTDIR=testing

SRCS=tool.cc support.cc
OBJS=$(subst .cc,.o,$(SRCS))

all: tool

tool: $(OBJS)
    $(CXX) $(LDFLAGS) -o tool $(OBJS) $(LDLIBS) 

tool.o: tool.cc support.hh

support.o: support.hh support.cc

clean:
    $(RM) $(OBJS)

distclean: clean
    $(RM) tool