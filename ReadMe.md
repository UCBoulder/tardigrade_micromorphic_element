######################
# tardigrade_micromorphic_element
######################
Developed by:

Nathan A. Miller (nami2227@colorado.edu)



Under the supervision of:

Dr. Richard Regueiro

Professor of Civil Engineering

The University of   
Colorado at Boulder



** Description **

An implementation of a micromorphic hex8 element. The form is 
intended to be used for easy implementation into a finite element 
or other variationally based framework. The code has been implemented 
into C++ code and utilizes the Eigen library. 

An implementation of a quasi-static finite element for use in the 
simulation code Abaqus is available in /src/cpp and can be built 
using:
> `make makefile_uel`

Note that this element uses the linear elastic material model but can 
use other material models by changing which is built in the makefile.
It is important to note that transient analyses in Abaqus are much 
more difficult to implement than was originally thought. This is due 
to restrictions on the degrees of freedom.

A second implementation which is more general and designed for use in 
the simulation code MOOSE can be built from code available in /src/cpp
 using the commands
> `make -f makefile_libmicromat`

> `make -f makefile_libmicrobalance`

This produces two shared librarys of micromorphic material models which 
can be linked to when compiling the MOOSE application "tardigrade" 
which is the implementation of micromorphic continuum mechanics.

** Upcoming work **

- Improve documentation
- Add computation of micro-gyration tensor
- Add plasticity models
- Add viscoelastic models

** Description of directories **

- .\doc: Where all of the genereated LaTeX documentation will be located. This includes 
       the main report of the testing of the code, the users manual, the programmers 
       manual and the theory manual. The resulting manuals are in PDF form.
       
- .\src: Where all of the source code for the micromorphic element lies. There are three
       subdirectories, python, cpp, and fortran which contain all python, C++, and 
       Fortran files required to compile the code.

** Python Requirements **

Users should compile the code within the MOOSE conda environment. If the python interface is to be used then
`cython` and `pytest` will also be needed. These can be installed via

    conda activate moose
    conda install cython pytest

** LaTeX Requirements **

The documentation requires an installation of LaTeX and Bibtex. It is assumed that the commands for these 
functions are pdflatex and bibtex respectively

*** LaTeX packages ***

- Report/Users Manual/Programmers Manual:
    - \usepackage{listings, xcolor, subcaption, placeins}
    - \usepackage{undertilde}
    - \usepackage{algorithm,algpseudocode}
    - \usepackage{multicol}
    - \usepackage{makecell}
    - \usepackage[table]{colortbl}

- Theory Manual:
    - beamer theme Pittsburg
    - \usepackage[utf8]{inputenc}
    - \usepackage{amsmath}
    - \usepackage{amsfonts}
    - \usepackage{amssymb}
    - \usepackage{undertilde}
    - \usepackage{bm}
    - \usepackage{subcaption}

** C++ Compiler Requirements **

Requires the library [Eigen](http://eigen.tuxfamily.org) which is a collection of header files and does
not require any compilation. The user must define the path to this library in `config.py`
(`eigen_location = /absolute/path/to/eigen`)

** Code Setup **

The code can be tested and the documentation generated by issuing the command >python setup.py
where "python" is the call to the python 2.7 installation. The results of the tests can be viewed in
the directory `.\doc\Report\`. Prior to running setup.py the relevant configuration questions 
in config.py.default should be set to the values for the system in question. Then a copy of 
config.py.default should be make in the same directory called config.py. Running >`python setup.py` 
will now configure the makefiles with the required definitions.

A full test of the code can be executed with the command
> `python run_tests.py -v`

where "python" is the call to the python 2.7 installation. This will run the unittest module 
which will execute the test functions in the code. The -v option allows more information to 
be printed to the screen. This is useful to understand which tests are being executed.

The code is also intended to generate the documentation. Currently, this is prevented because 
the LaTeX file generation is not particularly safe and can create errors. This will be 
investigated more soon.

** CPP Code **

Requires the GCC compiler (or other) though it defaults to gcc.

A local installation of GCC can be used by:

1. Downloading the gcc compiler and untaring it
2. Changing to the directory and running: `./contrib/download_prerequisites`
3. Running the command `./configure --prefix=/absolute/path/to/install/directory`
4. Running the command `make`
5. Before running a program compiled with this compiler one must set the environment 
   variable `export LD_LIBRARY_PATH=/absolute/path/to/install/directory`

Also requires the library [Eigen](http://eigen.tuxfamily.org) which requires that 
the path is defined in `config.py` (`eigen_location = /absolute/path/to/eigen`)

You will need to copy config.py.default to config.py. This allows git to track changes which 
are not system specific to be tracked in the default file while allowing the user to make 
changes locally. This same approach is also used in the makefiles with makefile.default being 
the file tracked by git.

** Abaqus environment variable setup **

The environment variable folder must be modified so that the user element library can be read 
and used in an analysis.
