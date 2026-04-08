import sys
import glob
from setuptools import setup, find_packages, Extension

c_sources = (
    glob.glob('NumbaQuadpack/src/*.c') +
    ['NumbaQuadpack/libcquadpack_module.c']
)

quadpack_ext = Extension(
    'NumbaQuadpack.libcquadpack',
    sources=c_sources,
    include_dirs=['NumbaQuadpack/include'],
    libraries=['m'] if sys.platform.startswith('linux') else [],
)

setup(
    name='ulysses',
    version='2.0.2',
    description='ULYSSES: Universal LeptogeneSiS Equation Solver',
    url='https://github.com/earlyuniverse/ulysses',
    author='Alessandro Granelli, Christopher Leslie, Yuber Perez Gonzalez, Brian Shuve, Holger Schulz, Jessica Turner, Rosie Walker',
    author_email='jessicaturner.5390@gmail.com',

    packages=find_packages(),
    include_package_data=True,

    install_requires=[
        'numpy>=1.22,<2.0',      # numba requires numpy<2.0
        'scipy>=1.9',
        'matplotlib>=3.5',
        'progressbar',
        'pandas',
        'python-dateutil',
        'termcolor',
        'tqdm',
        'numba>=0.56,<0.61',     # 0.56+ for cfunc stability; <0.61 until numpy2 support confirmed
        'mpmath',
        'mpltern',
        'multiprocess',
        'ipykernel'
    ],
    extras_require={
        'nest': ['pymultinest'],  # optional: not available on Windows
    },

    python_requires='>=3.9,<3.13',  # numba supports 3.9-3.12 as of 0.59
    scripts=['bin/uls-calc', 'bin/uls-scan', 'bin/uls-nest', 'bin/uls-models', 'bin/uls-scan2D'],

    ext_modules=[quadpack_ext],
)
