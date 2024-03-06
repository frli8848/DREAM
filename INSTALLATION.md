# Introduction

The DREAM Toolbox can be installed both using pre-compiled binaries and from source code. Binaries are
currently available for Linux (`x86_64`), Windows (`x86_64`), and for macOS (Intel `x86_64` Macs). The binaries
are compiled using generic compiler flags and should, therefore, run on most setups. If you want higher performance
then it is recommended that you compile DREAM from source as
described below.

# Binary Installation

There is (experimental) binaries of the MATLAB mex-files for Windows, Linux, and macOS here:

https://github.com/frli8848/DREAM/releases


The zip-files can be found under the _Assets_ link for each release.

## Windows MATLAB (mex) Binaries

Download the zip-file for Windows, uncompress it and move the files to a suitable location.
Also install the [MSVC++ redistributable][https://learn.microsoft.com/en-US/cpp/windows/latest-supported-vc-redist]
Then, follow the instructions on the _Post Installation Setup_ section below.

## Linux MATLAB Binaries

Download load the zip-file for Linux, uncompress it and move the files to a suitable location.
Then, follow the instructions on the _Post Installation Setup_ section below.

## macOS MATLAB Binaries

We are using Miniconda packages when we build the macOS binaries so one first needs to install Miniconda from here:

https://conda.io/projects/conda/en/latest/user-guide/install/macos.html

and then start a new Terminal (or Terminal tab) and install the FFTW lib:

```
$ sudo conda install -c conda-forge fftw
```

For macOS the automatic GitHub CI/CD builds puts `conda` packages in `/usr/local/miniconda` where as the
Miniconda installer puts the `conda` packages in `/opt/miniconda3`. And, one will get a "Library not loaded"
error if one uses the GitHub builds on a standard Miniconda install. One can fix this with a simple symlink:
```
 $ sudo ln -s /opt/miniconda3 /usr/local/miniconda
```

Then, follow the instructions on the _Post Installation Setup_ section below.

# Building the DREAM Toolbox Source Code

The  https://github.com/frli8848/DREAM.git repository contains the current development sources for the DREAM Toolbox.
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
to install a compiler tool chain and CMake for the corresponding platform.

The MATLAB builds are enabled with `-DBUILD_MEX=on` CMake flag (defaults to `off`)
and the Octave builds with`-DBUILD_OCT=on` flag (defaults to `on`).

## Linux

We have currently build the DREAM Toolbox on two Linux distributions: Ubuntu 22.04 LTS
and Gentoo Linux. We assume that a compiler tool chain is already installed, such as, `gcc`
or `clang`.

First install the FFTW library (not needed if Octave is installed) which can be done using:
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
DREAM/build $ cmake -DCMAKE_CXX_FLAGS="-O3 -march=native" -DBUILD_OCT=on -DBUILD_MEX=on ..
DREAM/build $ make -j8
```
which will build DREAM with both MATLAB and Octave support.

## macOS

On macOS, first install the macOS developer tools `Xcode` with the command line tools.
This will install  the `clang` compilers, `git` and other tools. Next install CMake from here: https://cmake.org/install/
and enable the command line tools from the Terminal:
```
$ sudo "/Applications/CMake.app/Contents/bin/cmake-gui" --install
Linked: '/usr/local/bin/cmake' -> '/Applications/CMake.app/Contents/bin/cmake'
Linked: '/usr/local/bin/ctest' -> '/Applications/CMake.app/Contents/bin/ctest'
Linked: '/usr/local/bin/cpack' -> '/Applications/CMake.app/Contents/bin/cpack'
Linked: '/usr/local/bin/cmake-gui' -> '/Applications/CMake.app/Contents/bin/cmake-gui'
Linked: '/usr/local/bin/ccmake' -> '/Applications/CMake.app/Contents/bin/ccmake'

```
Now, install FFTW using the Miniconda package manager (needed for MATLAB only builds)
```
 $ sudo conda install -c conda-forge fftw

```

And finally, build the optimized binaries for macOS by following the `git clone ...` and CMake build instructions
in the Linux build section above.

## Windows (mex Files)

First install the MSVC Toolchain (the community version works fine) and we assume that MATLAB is
already installed. Then install `git`, `cmake` and `Miniconda` or `Anaconda`:

* https://git-scm.com/download/win
* https://cmake.org/download/
* https://docs.conda.io/en/latest/miniconda.html or https://www.anaconda.com/download

When you install `cmake` add it to the paths for all users. Then start the "Anaconda Powershell prompt" and
run:
```
 >  conda install fftw
```

Now start a "Git Bash" shell and clone the DREAM sources using:
```bash
MINGW ~ $ git clone https://github.com/frli8848/DREAM.git
MINGW ~ $ cd DREAM
MINGW ~/DREAM (master) $ mkdir build
```

Then start an Anaconda Powershell and run
```
 > cd DREAM/build
 > cmake -DCMAKE_CXX_FLAGS="-O2 -EHsc" -DBUILD_MEX=on -DBUILD_OCT=off ..
 > cmake --build . --config Release
```
The `"-EHsc"` flag is for silencing MSVC exception handling warnings.

If everything builds then the newly build mex-files should be located in the folder:
```
 DREAM/build/Release
```

## Python (on Linux) - Experimental

First, one needs to install `pybind11`, `numpy` and `matplotlib` (to run the tests). On Gentoo Linux
```
# emerge dev-python/pybind11
# emerge dev-python/numpy
# emerge dev-python/matplotlib
```
and on Ubuntu Linux
```
$ sudo apt install python3-pybind11
$ sudo apt install python3-numpy
$ sudo apt install python3-matplotlib
```
Then (clone if needed) configure and build using (`mex` and `oct` builds are optional and can be switched off):
```bash
$ git clone https://github.com/frli8848/DREAM.git
$ cd DREAM
DREAM $ mkdir build && cd build
DREAM/build $ cmake -DCMAKE_CXX_FLAGS="-O3 -march=native" -DBUILD_OCT=on -DBUILD_MEX=on -DBUILD_PYTHON=on ..
DREAM/build $ make -j8
```

Finally, set the `PYTHONPATH` to your `build/python` folder (see the corresponding section below) and run the tests
using, for  example:
```
$ cd DREAM/python/tests
$ DREAM/python/tests $ python3 test_dreamrect.py
```
A window should appear with the SIR plots.

## Julia (on Linux) - Experimental

First install Julia, where on Ubuntu (22.04 LTS) there is no package so one have to do something like
```bash
$ wget https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.5-linux-x86_64.tar.gz
$ tar xf julia-1.8.5-linux-x86_64.tar.gz
$ sudo cp -a julia-1.8.5 /opt/
$ sudo ln -s /opt/julia-1.8.5/bin/julia /usr/local/bin/julia
```

On Gentoo Julia is in unstable so one have to add
```
sci-mathematics/dsfmt
sci-libs/openlibm
dev-lang/julia
app-emacs/julia-mode
```
to a file in `/etc/portage/package.accept_keywords/` and add the USE flag
```
net-misc/curl ssh
```
to a file in  `/etc/portage/package.use/`. Then one can install Julia using
```bash
$ sudo emerge dev-lang/julia
```

Now install the  `CxxWrap` package (and some other useful packages). Start Julia and press `]`
to enter `Pkg` mode

```bash
$ julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.8.5 (2023-01-08)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/

julia>
(@v1.8) pkg> add Plots
(@v1.8) pkg> add TickTock
(@v1.8) pkg> add DSP
(@v1.8) pkg> add https://github.com/barche/libcxxwrap_julia_jll.jl.git
(@v1.8) pkg> add CxxWrap
```

Then (clone if needed) configure and build using:
```bash
$ git clone https://github.com/frli8848/DREAM.git
$ cd DREAM
DREAM $ mkdir build && cd build
DREAM/build $ cmake -DCMAKE_CXX_FLAGS="-O3 -march=native" -DBUILD_JULIA=on ..
DREAM/build $ make -j8
```

# Building with GPU Acceleration Support - OpenCL Setup

We have experimental OpenCL (GPU) support for a few functions which currently includes:
* `das`
* `dreamrect`
* `dreamcirc`
* `rect_sir`

The support is enabled by using the cmake option `-DUSE_OPENCL=on`. One needs driver support for
the type of GPU used (Nvidia, AMD, or Intel) which is normally provided by the GPU vendor tools
like, for example, CUDA for Nvidia or the ROCm stack from AMD.

## Linux

On Gentoo first install the vendor specific OpenCL tools, for example CUDA (or ROCm etc.),
the OpenCL C++ headers, and and then the `xxd` tool from vim (which we use to convert OpenCL
`.cl` sources to header files):
```bash
# emerge dev-util/nvidia-cuda-toolkit
# emerge dev-util/opencl-headers dev-libs/clhpp
# emerge app-editors/vim-core

```
or on Ubuntu,
```bash
# sudo apt install nvidia-cuda-toolkit-gcc
# sudo apt install opencl-c-headers opencl-clhpp-header
```

Now setup the build as described in the Linux build section above and add the cmake option `-DUSE_OPENCL=on`
when configuring the build like, for example,
```bash
-snip-

DREAM/build $ cmake -DCMAKE_CXX_FLAGS="-O3 -march=native" -DBUILD_OCT=on -DBUILD_MEX=on -DUSE_OPENCL=on ..

-snip-
```
Then follow the build instructions in the Linux build section above.

## macOS

Nothing here yet.

## Windows

* Install CUDA (Nvidia), ROCm (AMD) [not tested], or Intel tools [not tested] depending of which GPU is in use.

* Install the OpenCL headers. Start `Git Bash` and run
  ```
  $ git clone https://github.com/KhronosGroup/OpenCL-Headers.git
  $ git clone https://github.com/KhronosGroup/OpenCL-CLHPP.git
  ```
  Then create the folder `CL` in the `DREAM/include` folder and copy all OpenCL header
  files (both `.h` and `.hpp`) to that folder; they are located in `OpenCL-Headers/CL`
  and `OpenCL-CLHPP/include/CL`, respectively.

* Start an Anaconda powershell and install the `xxd` tool (part of vim):
  ```
  > conda install --channel=conda-forge vim
  ```
* Finally, using an Anaconda Powershell, configure and build using:
  ```
 > cmake -DCMAKE_CXX_FLAGS="-O2 -EHsc" -DBUILD_MEX=on -DBUILD_OCT=off -USE_OPENCL=on ..
 > cmake --build . --config Release
 ```

# Post Installation Setup

## MATLAB and Octave

After building/installing the MATLAB mex-files and/or Octave oct-files add the build folder, or biniary install folder, to
the MATLAB/Octave path. That is, for Octave
add
```
addpath('/<YOUR-HOME-DIR>/<PATH-TO-DREAM-SOURCES>/DREAM/build')
```
to the `~/.octaverc` file, and for MATLAB add it to the `~/Documents/MATLAB/startup.m` file.

## Python

After building the Python bindings add the build/installation folder to your `PYTHONPATH`
using, for example
```bash
$ export PYTHONPATH="$PYTHONPATH:$HOME/DREAM/build/python"
```
if you have your DREAM sources in `$HOME/DREAM`.

## Julia

After building the Julia bindings add the build/installation folder to the Julia `LOAD_PATH`
by adding it to your  `~/.julia/config/startup.jl` by adding the line:
```
push!(LOAD_PATH, "/<YOUR-HOME-DIR>/<PATH-TO-DREAM-SOURCES>/DREAM/build/julia")
```

You should now be able to run, for example
```bash
DREAM/julia/tests $ julia -i test_dreamrect.jl
```
