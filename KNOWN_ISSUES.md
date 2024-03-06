# MATLAB

## Linux

* To avoid (work-around) a conflict with the FFTW lib used by DREAM and the one bundled with MATLAB one can
(on Linux systems) preload the one used by DREAM using, for example
```bash
 $ export LD_PRELOAD=/usr/lib64/libfftw3.so
```
The path to `libfftw3.so` may differ between different Linux distributions

* When OpenCL support is enable one can get an error similar to:
```
/usr/opt/MATLAB/R2020b/bin/glnxa64/../../sys/os/glnxa64/libstdc++.so.6: version `GLIBCXX_3.4.29' not found
(required by /home/fl/DREAM/build/das.mexa64
```
which is due to a conflict with `libstdc++.so.6`  bundled with MATLAB and the one used by the Linux system.
A work-around is to `LD_PRELOAD` the system `libstdc++.so.6` using, for example.
```bash
$ export LD_PRELOAD=$LD_PRELOAD:/usr/lib/gcc/x86_64-pc-linux-gnu/12/libstdc++.so.6
```
The path to `libstdc++.so` may differ between different Linux distributions and between different compilers.

## CTRL-C Support

Even though we have code to catch CTRL-C interrupts (to force a stop of a lengthy simulation) it is not used
since the MATLAB, Octave, and Python interpreters also catch the CTRL-C signal (before the DREAM code can register
the CTRL-C handler). On Julia we can, however, switch that off using:
```
# This disables Julia's (REPL) signal handler for SIGINT
# so that we can CTRL-C and interrupt our code if we want.
ccall(:jl_exit_on_sigint, Cvoid, (Cint,), 0)

```
