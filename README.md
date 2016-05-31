# README #

### TO DO LIST PHYSICS ###

1. Add self gravity
1. Compute Shear Heating and add it to temperature
1. Compute Viscous flow law according to the viscous strain rate, as opposed to the total strain rate
1. Add compressibility

### TO DO LIST OPTIMIZATION ###

1. optimize marker to cells to make less boolean tests
1. add subgrid diffusion on markers
1. Physics_InterpFromParticlesToCell: values of T might need to be overwritten for periodic nodes
1. Create Numbering->boundary and Numbering->boundaryNeighbour to optimize the looping the boundaries for EC arrays
1. For periodic boundaries instead of numbering lines like this: 0 1 2 3 0 1, change to numebring: 3 0 1 2 3 1, or when computing interpolation to Boundaries use the solution from the inner node
1. Remove particles from the air
1. implement swiss cross grid
1. inject particle when cell receives no contribution

### TO DO LIST VISUALIZATION ###
hCalling Python to do data treatment would be ideal. See :ttp://www.linuxjournal.com/article/8497
See also:
https://flamingoengine.wordpress.com/2010/06/28/pyopengl-too-slow-use-c/

### TO DO LIST SETUP ###
1. Add sandbox type BC including constant heat flux boundaries
1. Read geomIO file and perform point in polygon test

### DONE ###
- Stokes solver
- Non linear rheology
- Penalty method
- Pure shear BC
- Periodic BC
- Linked list of markers
- Velocity, strain rate and viscosity visualization
- Interpolation of Viscosity from markers to cell centers
- Interpolation of Viscosity froms cell centers to markers
- Interpolation routine from cell centers to nodes
- Add gravity
- Use texture instead of triangles
- Visualize markers using a geometry shader or point sprites
- Add passive grid
- Add temperature
- Add elasticity
- Add plasticity


## INSTALLATION

# Required libraries
Pardiso
omp

# Optional library for graphics

To use opengl:

- glfw
- glew

to save images:

- libgpng

external but useful package:

- ffmpeg () to create videos

## Installation on MacOSX

# Using Homebrew package manager
- brew install gcc-5 // (by default MacOS calls clang which doesn't work ery well with omp)
- brew install glfw3
- brew install glew
- brew install libpng
- brew install ffmpeg --with-fdk-aac --with-ffplay --with-freetype --with-libass --with-libquvi --with-libvorbis --with-libvpx --with-opus --with-x265


If you get this error message at execution:

dyld: Library not loaded: /usr/local/lib/libgfortran.3.dylib
  Referenced from: /usr/local/lib/libpardiso500-MACOS-X86-64.dylib
  Reason: image not found

then, manually copy the following libraries from /usr/local/lib/gcc/5/ to /usr/local/lib: 

- libfortran.3.dylib
- libgfortran.dylib
- libgfortran.a
- libgfortran.spec
- libgomp.1.dylib
- libgomp.dylib
- libgomp.a
- libgomp.spec
- libquadmath.0.dylib
- libquadmath.dylib
- libqadmath.a


To run the code with eclipse waking into account the environment variable OMP_NUM_THREADS (required by pardiso), follow this: http://stackoverflow.com/questions/829749/launch-mac-eclipse-with-environment-variables-set


## Tips

Encoding a video readable with quicktime using ffmpeg:

ffmpeg -i Frame_%05d.png -f mp4  -vcodec h264 -pix_fmt yuv420p  Movie.mp4


Note to install matplotlib and pyopengl

pygame:
install python 3

follow this:
http://brysonpayne.com/2015/01/10/setting-up-pygame-on-a-mac/

install matplotlib: 

follow : http://matplotlib.org/1.3.1/users/installing.html

but in the last step use: sudo python3 setup.py install