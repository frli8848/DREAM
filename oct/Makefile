#!/usr/local/bin/make

# Copyright (C) 2006,2007,2008 Fredrik Lingvall
#
# This file is part of the DREAM Toolbox.
#
# The DREAM Toolbox is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2, or (at your option) any
# later version.
#
# The DREAM Toolbox is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU General Public License
# along with the DREAM Toolbox; see the file COPYING.  If not, write to the
# Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.

# $Revision: 826 $ $Date: 2015-04-29 08:05:47 +0200 (Wed, 29 Apr 2015) $ $LastChangedBy: frli8848 $

###################################
#
# Makefile for building DLD Octave
# files from the C-sources.
#
###################################

include ../Make.inc

#
# Compiler for the attenuation code (see also att_pthread.o below).
#

# Needs a C99 compiler (gcc for example).
ATT = att.o

# Non iso C99 compiler (uses C++) and should work with MSVC.
#ATT = att_octave.o

##################################################################
#
# Note:
#
# Octave must be compiled with (at least):
#
#	export CXXFLAGS=-pthread (export CFLAGS==-pthread)
#
# and
#	./configure  --enable-shared --enable-dl
#
##################################################################

DLDCC		= $(CXX)
#DLDCC		= mkoctfile

# Octave is compiled with -fPIC and when using -shared the DLD (.oct) functions
# must therefore also be compiled with -fPIC (according to the gcc manual).
FPIC		= -fPIC


DREAM_CFLAGS	= $(CFLAGS)  $(ATT_FFTW)
DREAM_CXXFLAGS	= $(CXXFLAGS)  $(ATT_FFTW) -std=c++11


# Octave options for g++/gcc linking.
DLD_FLAGS	= -shared -Wl,-Bsymbolic
#DLD_FLAGS	=

#
# LD flags (from gcc manual).
#

#        Note - if the linker is being invoked indirectly, via a compiler driver
#        (eg gcc) then all the linker command line options should be prefixed by
#        -Wl, (or whatever is appropriate for the  particular  compiler  driver)
#        like this:
#
#                  gcc -Wl,--startgroup foo.o bar.o -Wl,--endgroup
#
#        This  is  important,  because otherwise the compiler driver program may
#        silently drop the linker options, resulting in a bad link.
#
#       -Bsymbolic
#            When  creating  a shared library, bind references to global symbols
#            to the definition within the shared library, if any.  Normally,  it
#            is  possible for a program linked against a shared library to over-
#            ride the definition within the shared library.  This option is only
#            meaningful on ELF platforms which support shared libraries.
#

DREAM_INCLUDE		= -I../include $(OCTAVE_HEADERS)

# Needed for *_p functions (you may need to build Octave with -pthread flag). TODO: is this needed with C++11?
PARLIBS		= -lpthread

LIBDIRS		= $(OCT_LIBS) $(FFW_LIB)

#
# Rules for the obj. files.
#

.c.o:
	$(CC) -DOCTAVE $(DREAM_CFLAGS) $(FPIC) $(DREAM_INCLUDE) -c $<
#	$(DLDCC) -DOCTAVE $(DREAM_CFLAGS) $(FPIC) $(DREAM_INCLUDE) -c $<


.cc.o:
	$(CXX) -DOCTAVE $(DREAM_CXXFLAGS) $(FPIC) $(DREAM_INCLUDE) -c $<
#	$(DLDCC) -DOCTAVE $(DREAM_CXXFLAGS) $(FPIC) $(DREAM_INCLUDE) -c $<

all: \
	dreamline.oct \
	dreamrect.oct \
	dreamrect_f.oct \
	dreamcirc.oct \
	dreamcirc_f.oct \
	dreamsphere_f.oct \
	dreamsphere_d.oct \
	dreamcylind_f.oct \
	dreamcylind_d.oct \
	dream_arr_rect.oct \
	dream_arr_circ.oct \
	dream_arr_cylind_f.oct \
	dream_arr_cylind_d.oct \
	dream_arr_annu.oct \
	dream_apodwin.oct \
	dream_att.oct

#
# Compile att.c with pthread support (uses pthread mutexes).
#

# C99
att_pthread.o:
	$(CC) $(DREAM_CXXFLAGS) $(FPIC) $(DREAM_INCLUDE) -DPARALLEL -c att.c -o $@
#	$(DLDCC) $(DREAM_CXXFLAGS) $(FPIC) $(DREAM_INCLUDE) -DPARALLEL -c att.c -o $@


# C++
#att_pthread.o:
#	$(CXX) $(DREAM_CXXFLAGS) $(FPIC) $(DREAM_INCLUDE) -DPARALLEL -c att_octave.cc -o $@


#
# Rules for the OCT-files.
#

# Single element transducers

dreamline.oct : oct_dreamline.o dreamline.o att_pthread.o fft_octave.o dream_error.o
	$(DLDCC) $(DLD_FLAGS) $(PARLIBS) $(LIBDIRS)  $^ -o $@

dreamrect.oct : oct_dreamrect.o dreamrect.o  att_pthread.o fft_octave.o dream_error.o
	$(DLDCC) $(DLD_FLAGS) $(PARLIBS) $(LIBDIRS)  $^ -o $@

dreamrect_f.oct : oct_dreamrect_f.o dreamrect_f.o  att_pthread.o fft_octave.o dream_error.o
	$(DLDCC) $(DLD_FLAGS) $(PARLIBS) $(LIBDIRS)  $^ -o $@

dreamcirc.oct : oct_dreamcirc.o dreamcirc.o  att_pthread.o fft_octave.o dream_error.o
	$(DLDCC) $(DLD_FLAGS) $(PARLIBS) $(LIBDIRS)  $^ -o $@

dreamcirc_f.oct : oct_dreamcirc_f.o dreamcirc_f.o  att_pthread.o fft_octave.o dream_error.o
	$(DLDCC) $(DLD_FLAGS) $(PARLIBS) $(LIBDIRS)  $^ -o $@

dreamsphere_f.oct : oct_dreamsphere_f.o dreamsphere_f.o  att_pthread.o fft_octave.o dream_error.o
	$(DLDCC) $(DLD_FLAGS)  $(PARLIBS) $(LIBDIRS)  $^ -o $@

dreamsphere_d.oct : oct_dreamsphere_d.o dreamsphere_d.o  att_pthread.o fft_octave.o dream_error.o
	$(DLDCC) $(DLD_FLAGS) $(PARLIBS) $(LIBDIRS)  $^ -o $@

dreamcylind_f.oct : oct_dreamcylind_f.o dreamcylind_f.o  att_pthread.o fft_octave.o dream_error.o
	$(DLDCC) $(DLD_FLAGS) $(PARLIBS) $(LIBDIRS)  $^ -o $@

dreamcylind_d.oct : oct_dreamcylind_d.o dreamcylind_d.o  att_pthread.o fft_octave.o dream_error.o
	$(DLDCC) $(DLD_FLAGS) $(PARLIBS) $(LIBDIRS)  $^ -o $@

# Arrays

dream_arr_rect.oct : oct_dream_arr_rect.o dream_arr_rect.o arr_functions.o att_pthread.o fft_octave.o dream_error.o
	$(DLDCC) $(DLD_FLAGS) $(PARLIBS) $(LIBDIRS)  $^ -o $@

dream_arr_circ.oct : oct_dream_arr_circ.o dream_arr_circ.o arr_functions.o att_pthread.o fft_octave.o dream_error.o
	$(DLDCC) $(DLD_FLAGS) $(PARLIBS) $(LIBDIRS)  $^ -o $@

dream_arr_cylind_f.oct : oct_dream_arr_cylind_f.o dream_arr_cylind_f.o arr_functions.o att_pthread.o fft_octave.o dream_error.o
	$(DLDCC) $(DLD_FLAGS) $(PARLIBS) $(LIBDIRS)  $^ -o $@

dream_arr_cylind_d.oct : oct_dream_arr_cylind_d.o dream_arr_cylind_d.o arr_functions.o att_pthread.o fft_octave.o dream_error.o
	$(DLDCC) $(DLD_FLAGS) $(PARLIBS) $(LIBDIRS)  $^ -o $@

dream_arr_annu.oct : oct_dream_arr_annu.o dream_arr_annu.o att_pthread.o fft_octave.o dream_error.o
	$(DLDCC) $(DLD_FLAGS) $(PARLIBS) $(LIBDIRS)  $^ -o $@

#
# Misc routines.
#

dream_apodwin.oct : oct_dream_apodwin.o arr_functions.o dream_error.o
	$(DLDCC) $(DLD_FLAGS) $(LIBDIRS) $^ -o $@

dream_att.oct : oct_dream_att.o $(ATT) fft_octave.o dream_error.o
	$(DLDCC) $(DLD_FLAGS) $(LIBDIRS) $^ -o $@



install:
	cp -a *.oct ../$(INSTALL_DIR)

.PHONY : clean

clean:
	rm -f  *.o *~ *.obj *.rsp

distclean:
	rm -f   *.o *~ *.obj *.rsp *.oct
