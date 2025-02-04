from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy
from numpy import get_include

ext = Extension("cython_process", sources=["cython_process.pyx"], include_dirs=['.', get_include()])
setup(name="cython_process", ext_modules=cythonize([ext]),include_dirs=[numpy.get_include()])
