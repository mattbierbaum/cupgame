import  os
from setuptools import setup
from distutils.extension import Extension
from distutils.command.build_ext import build_ext

ext = Extension('_cupgamelib', sources=['src/cupgamelib.i', 'src/cupgamelib.c'],
                libraries=['m'], extra_compile_args=["-O2", "-funroll-loops", "-Wall"],
                include_dirs=['src'])

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
