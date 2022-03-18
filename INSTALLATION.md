# Introduction

The DREAM Toolbox can be installed both using pre-compiled binaries and from source code. Binaries are
currently available for Linux (x86_64) and for macOS (Intel x86_64 Macs). The binaries are compiled
using generic compiler flags and should, therefore, run on most setups. If you want higher performance
then it is recommended that you compile DREAM from source as
described below

# Binary Installation

There is (experimental) binaries of the Matlab mex-files for Linux and macOS here:

https://github.com/frli8848/DREAM/releases

## Linux Matlab Binaries

Nothing here yet.

## macOS Matlab Binaries

We are using Miniconda pcakages when we build the macOS binaries so one first needs to install Miniconda from here:

https://conda.io/projects/conda/en/latest/user-guide/install/macos.html

and then start anew Terminal (or Terminal tab) and install the FFTW lib:

```
$ sudo conda install -c conda-forge fftw
```

For macOS the automatic GitHub CI/CD builds puts `conda` packages in `/usr/local/miniconda` where as the
Miniconda installer puts the `conda` packages in `/opt/miniconda3`. And, one will get a "Library not loaded"
error if one uses the GitHub builds on a standard Miniconda install. One can fix this with a simple symlink:
```
 $ sudo ln -s /opt/miniconda3 /usr/local/miniconda
```

## Windows Matlab Binaries

Nothing here yet.

# Building the DREAM Toolbox Source Code

The  https://github.com/frli8848/DREAM.git repsitory contains the current development sources for the DREAM Toolbox.
One can obtain the source code using:
```
 $ git clone https://github.com/frli8848/DREAM.git
```
or
```
 $ git clone git@github.com:frli8848/DREAM.git
```

The old legacy (and now obsolete) code can be
found here: https://sourceforge.net/projects/dreamtoolbox/

We aim to build the code and documentation on all platforms using `cmake` and one needs
to install  a compiler tool chain and CMake for the corresponding platform.

The Matlab builds are enabled with `-DBUILD_MEX=on` CMake flag (defaults to `off`)
and the Octave builds with`-DBUILD_OCT=on` flag (defaults to `on`).

## Linux

We have currently build the DREAM Toolbox on two Linux distributions: Ubuntu 20.04 LTS
and Gentoo Linux. We assume that a complier toolchain is already installed, such as, `gcc`
ocr `clang`.

First install the FFTW library which can be done using:
```
 $ sudo apt -yq update
 $ sudo apt install libfftw3-dev
```
on Ubuntu Linux, or
```
$ sudo emerge sci-libs/fftw

```
on Gentoo Linux. To Install (a resent) Octave version one can do

```
$ sudo add-apt-repository ppa:devacom/science
$ sudo apt install octave
$ sudo apt install liboctave-dev
```
on Ubuntu, and on Gentoo Linux, add `sci-mathematics/octave` to a file in the
`/etc/portage/package.accept_keywords/` folder and install Octave using
```
# emerge sci-mathematics/octave

```
Installing Octave will also pull in FFTW so the FFTW install steps above is then not needed.

For other Linux distributions consult your package manger's documentation.

Now, to build optimized Linux binaries do, for example,

```bash
$ git clone https://github.com/frli8848/DREAM.git
$ cd DREAM
DREAM $ mkdir build && cd build
DREAM/build $ cmake -DCMAKE_C_FLAGS="-O3 -march=native"  -DCMAKE_CXX_FLAGS="-O3 -march=native" -DBUILD_OCT=on -DBUILD_MEX=on ..
DREAM/build $ make -j8
```
which will build DREAM with bothe Matlab and Octave support.


## macOS

on macOS using the Miniconda package manager.

```
 $ sudo conda install -c conda-forge fftw

```

To build optimized binaries for macOS follow the same `git clone ...` and CMake build instructions
in the Lnux build section above

## Windows

Nothing here yet.

# Post Installation Setup

After building/installing the Matlab mex-files and/or Octave oct-files add the build folder, or biniary install folder, to
the Matlab/Octave path. That is, for Octave
add
```
addpath('/<YOUR-HOME-DIR>/<PATH-TO-DREAM-SOURCES>/DREAM/build')
```
to the `~/.octaverc` file, and for Matlab add it to the `~/Documents/MATLAB/startup.m` file.
