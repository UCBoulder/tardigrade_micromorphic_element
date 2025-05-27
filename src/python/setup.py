import os
from distutils.core import setup
from distutils.extension import Extension

import numpy
from Cython.Distutils import build_ext

project_name = 'micromorphic'
cxx_standard = '11'

###############################
# Get the include directories #
###############################

source_files = [os.path.abspath('../cpp/material_python_interface.cpp'),
                os.path.abspath('materials.pyx')]

include_dirs = [numpy.get_include(), os.path.abspath('../cpp/')]
library_dirs = [os.path.abspath('../../build/src/cpp')]

shared_libraries = ['micromat']

ext_modules = [Extension(project_name,
               sources=source_files,
               language='c++',
               libraries=shared_libraries,
               include_dirs=include_dirs,
               extra_compile_args=[f"-std=c++{cxx_standard}"],
               extra_link_args=[f"-L{d}" for d in library_dirs])]

for e in ext_modules:
    e.cython_directives = {'language_level': "3"}

setup(
  name = project_name,
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
