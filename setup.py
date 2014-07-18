# -*- coding: utf-8 -*-
'''This script compiles the cython that describes the kernel function, k_2d'''
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy


ext_modules= [ 
    Extension("k_2d", ["k_2d.pyx"],
              libraries=["m"],
              include_dirs = [numpy.get_include()]
              )
              ]

setup(
    name = 'k_2d',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,

    )
