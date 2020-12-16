# README #

Sunaba is a numerical code that can be used to model geodynamical thermomechanical processes. The code solves the Stokes and heat equations using a parallel implicit solver based on the finite-difference marker-in-cell method. It can simulate elasto-visco-plastic materials with non-linear temperature-dependent viscosity. The code is written in C for the calculation, and in Python for the user-interface. The code also provides real-time visualization of simulation results via OpenGL.
The details of the equations and algorithm are given in 

[Bauville, A., Furuichi, M., &Gerbault, M. I. (2020). Control of fault weakening on the structural styles of underthrusting-dominated non-cohesive accretionary wedges. Journal of Geophysical Research: SolidEarth, 125, e2019JB019220. https://doi.org/10.1029/2019JB019220](https://www.researchgate.net/publication/339624541_Control_of_Fault_Weakening_on_the_Structural_Styles_of_Underthrusting-Dominated_Non-Cohesive_Accretionary_Wedges)

## INSTALLATION
Compile the code from source using one of the makefile contained in the `Release` or `Debug` folder, or your own makefile.
To run the code you need an input file. e.g. to run the code from the `Release` folder using the input file in Setups: `./StokesFD ../Setups/input.json`
`input.json` files are generated using the python scripts. Examples of python scripts are given in the `Setups/Benchmarks` folder.

# Required software and libraries
- a C compiler (I recommend GCC or ICC rather than clang, because clang does not support openmp)
- the [Pardiso](https://www.pardiso-project.org/) solver
- Python and scientific computing libraries (e.g. [Anaconda](https://www.anaconda.com/))

# Optional library for graphics

For parallel calculation:

- omp

To use opengl:

- glfw
- glew

to save images:

- libgpng

external but useful package:

- ffmpeg () to create videos

# Installation tips for MacOSX
## Getting command line tools (including GCC) 
1. install Xcode from the AppStore
1. open Xcode, agree to the license
1. open a terminal
1. type $ xcode-select --install, then click install when prompted

note: the last step needs to be repeated in case of upgrading MacOSX
the command line tools include only GCC 4.2.1. For optimal speed, a newet version is recommended. It can be obtained from a package manager. 
Below, I use the package manager Homebrew

## Getting extra packages from Homebrew 
1. Download and install Homebrew on https://brew.sh
- brew install gcc // (by default MacOS calls clang which doesn't work very well with omp)
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


to restore the ownership of the folder to the file system:

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


To run the code with eclipse taking into account the environment variable OMP_NUM_THREADS (required by pardiso), follow this: http://stackoverflow.com/questions/829749/launch-mac-eclipse-with-environment-variables-set


# Other tips

Encoding a video readable with quicktime using ffmpeg:

ffmpeg -i Frame_%05d.png -f mp4  -vcodec h264 -pix_fmt yuv420p  Movie.mp4

to slow down a video: ffmpeg -i input.mkv -filter:v "setpts=0.5*PTS" output.mkv

if the image as an odd width or height. A error message will be outputted. A work around is to crop the images using imagemagick. 
To crop all the images in a folder to final Dimension W,H:

mogrify -crop WxH+0+0 *.png



batch cropping with image magick:

convert *.png -crop 300x500+100+70 result/cropped_image.png

Valgrind debug line:
valgrind --log-file="valgrindDebug.txt" --leak-check=full --show-leak-kinds=all ./StokesFD

-----

StokesFD includes custom tasks and shortcuts for Visual Studio Code.

To get the intellisense working you might need to edit "includePath" in the file c_cpp_properties.


