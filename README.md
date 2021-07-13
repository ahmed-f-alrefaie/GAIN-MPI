# GAIN-MPI

## Compiling on Wilkes2

To compile for use on Wilkes2, first make sure you are in the
gpu login nodes, you can do this by logging in like so:

```
ssh <diracusername>@login-gpu.hpc.cam.ac.uk
```

Clone the code and cd into the directory:

```
git clone https://github.com/ahmed-f-alrefaie/GAIN-MPI.git
cd GAIN-MPI
```

Then use the wilkes to makefile and make:

```
cp makefile.wilkes2 makefile
make -j10
```

You should now have the GAIN-MPI-KEPLER.x executable.
