# README #


### Holiday TO DO LIST ###

1. Redesign the injection so that it works good for local interpolation as well
1. Add sedimentation routines
1. Correct the Crank Nicholson routines to account for varying time step sizes
1. Check only the resiudals for the velocity


### Short term TO DO LIST ###

1. Debug the velocity advection
1. assign and check the pre-processor switches from the python input files. Recompile if necessary
1. Add
1. Add subgrid diffusion for phi and DeltaP0
1. make a better interface for boundary conditions and enforce 0 shear stress on particles near a free slip wall



### TO DO LIST PHYSICS ###

1. Add self gravity
1. Add compressibility

### TO DO LIST OPTIMIZATION ###

1. Debug the swiss cross grid
1. optimize marker to cells to make less boolean tests
1. Remove particles from the air


### TO DO LIST VISUALIZATION ###
Calling Python to do data treatment would be ideal. See :ttp://www.linuxjournal.com/article/8497
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
- Basic reading file and plotting in python
- Create a proper input system
- Clean the code: remove the dirty Darcy, clean the visualization
- Compute Viscous flow law according to the viscous strain rate, as opposed to the total strain rate
- add subgrid diffusion on markers
- Compute Shear Heating and add it to temperature
- Physics_InterpFromParticlesToCell: values of T might need to be overwritten for periodic nodes
- For periodic boundaries instead of numbering lines like this: 0 1 2 3 0 1, change to numebring: 3 0 1 2 3 1, or when computing interpolation to Boundaries use the solution from the inner node
- inject particle when cell receives no contribution
- Create a file saving system
- Add Darcy
- For periodic BC, contribution from the other side should be considered for Cells2Particles interpolation
- Find why Pe is not properly limited by Py
- Find the bug in the elasticity (soolid is compressible, so deviatoric strain rate should be used, really)
- Verify the advection of phi
- Put the colorscale as a user input
- Implement the proper deviatoric strain rate
- Finish implementing the phase list
- Implement better rheology
- Add sediments to the corner flow
- Strain softening
- Initial conditions for T and Darcy
- Give the same thermal diffusion to the air than to rocks, but put the temperature in the air to 0

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
- brew install gcc-5 // (by default MacOS calls clang which doesn't work very well with omp)
- brew install glfw3
- brew install glew
- brew install libpng
- brew install ffmpeg --with-fdk-aac --with-ffplay --with-freetype --with-libass --with-libquvi --with-libvorbis --with-libvpx --with-opus --with-x265

note: to upgrade homebrew on OSX10.12: 

if you get the error message:

/usr/local/Library/brew.sh: line 32: /usr/local/Library/ENV/scm/git: No such file or directory

then:

cd "$(brew --repository)" && git fetch && git reset --hard origin/master


if you get the error:

Error: /usr/local must be writable!

then use the following to switch the ownership of the folder (can be restored afterwards):

sudo chown -R $(whoami) /usr/local


to restor the ownership of the folder to the file system:

sudo chown root:wheel /usr/local



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

to slow down a video: ffmpeg -i input.mkv -filter:v "setpts=0.5*PTS" output.mkv



batch cropping with image magick:

convert *.png -crop 300x500+100+70 result/cropped_image.png


Nice python and Jupyter tutorial:
http://bebi103.caltech.edu/2015/tutorials.html

Valgrind debug line:
valgrind --log-file="valgrindDebug.txt" --leak-check=full --show-leak-kinds=all ./StokesFD