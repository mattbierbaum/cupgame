from setuptools import setup
from distutils.extension import Extension
import numpy as np

ext = Extension('_cupgamelib', sources=['src/cupgamelib.i', 'src/cupgamelib.c'],
                libraries=['m'], extra_compile_args=["-O2", "-funroll-loops", "-Wall", "-march=native"],
                include_dirs=['src', np.get_include()])

setup(
    name='cupgamelib',
    author='Matt Bierbaum',
    version='0.1',
    license="MIT",

    # this is necessary so that the swigged python file gets picked up
    py_modules=['cupgamelib', 'cupgame'],
    package_dir={'': 'src'},

    ext_modules = [ext],

    # since the package has c code, the egg cannot be zipped
    zip_safe=False
)
