from distutils.core import setup, Extension
from distutils.sysconfig import get_python_inc, get_python_lib
import os
import sys
import numpy

###################################################################
# build the extension
#

define_macros = []
undef_macros = []
extra_compile_args = []
include_dirs = []

libraries = []
library_dirs = []

# Use NumPy instead of Numeric or numarray
make_extension = Extension
include_dirs.append(numpy.get_include())
undef_macros.append('USE_NUMARRAY')

if not os.name == "posix":
    raise Exception, "os not supported"

ext = make_extension('_eventdft_c',
                     ['eventdft_wrap.c', 'period.c', 'utils.c'],
                     include_dirs=include_dirs,
                     libraries=libraries,
                     library_dirs=library_dirs,
                     define_macros=define_macros,
                     extra_compile_args=extra_compile_args)

###################################################################
# the package
#

setup(name="eventdft",
      version="0.99",
      description="Python interfaces to event_DFT",
      author="Scott Ransom",
      author_email="sransom@nrao.edu",
      packages=[''],
      package_dir={'eventdft':'eventdft_src'},
      ext_modules=[ext])
