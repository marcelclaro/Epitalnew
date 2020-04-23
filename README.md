Epitalnew
=======


This is a tentative to organize and make my code to calculate quantum levels, wave functions and some properties of semiconductor heterostructures (mainly grown by epitaxy) user-friendly and open source. Note that at beginning

Initially it was focused on  III-V materials, expanded to II-VI. Some VdW materials support is planned.

Since share it was not planned from the beginning, at some point I mislead some documentation. It will be updated as soon as possible.
If you want to use it extensively, please consider contribute to the documentation and with the code. If you need help send me a message.


License
-----------


![Image](https://licensebuttons.net/l/by-nc-sa/3.0/88x31.png "license")
This code is published with creative Commons license CC BY-NC-SA. It means This license lets others remix, adapt, and build upon your work non-commercially, as long as they credit and license their new creations under the identical terms.

Please, give me the acknowledgment if you publish something which use this code.


Requirements
-----------

I recommend to use CMake to compile it (if you are using VTK it is almost mandatory).

It uses gnuplot to plot some results.

It will be nice to have OpenMP since it has multi-thread support through it, and most of the operations scale with the number of processors.

It requires [armadillo](http://arma.sourceforge.net/) and [Eigen3](http://eigen.tuxfamily.org/) for some linear algebra and matrix operations.

HDF5 to save structured data.

[VTK](https://vtk.org/) is required for load some 3D models (Heterostructure3D.hpp). Sometimes is trick to install it. If you have problems and don't want to use this function, quit Heterostructure3D.hpp/Heterostructure3D.cpp files

How to use
-----------

The folder /src/Examples has some examples of tested simulations.

Check: /src/Examples/transfermatrix_example.cpp
It has an example how to calculate the energy and wavefunctions of an AlGaAs/GaAs quantum well.

To compile:

1. Create a folder inside /src with a name .build
2. Add your code
3. Include your main file on the CMakeLists.txt code
4. Run CMake and check if all the dependencies are satisfied
5. Run "make epital" and Good Luck!

Example:

```
cd src
(write our code and include it on the CMakeLists)
mk .build
cd .build
cmake ..
make
```


Dependencies to install on Ubuntu 18.04:
```
sudo apt install cmake
sudo apt install gnuplot
sudo apt install libarmadillo-dev
sudo apt install hdf5-tools
sudo apt install libeigen3-dev
sudo apt install vtk7
```
If something still missing, try [auto-apt](http://manpages.ubuntu.com/manpages/trusty/man1/auto-apt.1.html)

If you use Windows, start to think in use Linux.


TODO
-----------
* Update Docs
* Do test routines
* fix the files on Fixme
* Improve examples
* VdW materials InSe, GaSe and MoSe2
