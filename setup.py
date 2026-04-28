import os
import shutil
import sys
import glob
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext as _build_ext


class build_ext(_build_ext):
    """Copy the compiled extension into the source tree after building.

    This ensures driver.py (which looks next to itself) finds the shared
    library regardless of install mode — regular, editable, or in-place —
    even when Python is run from inside the repo root.
    """
    def run(self):
        super().run()
        for ext in self.extensions:
            built = self.get_ext_fullpath(ext.name)
            if not os.path.exists(built):
                continue
            dest = os.path.join('NumbaQuadpack', os.path.basename(built))
            if os.path.abspath(built) != os.path.abspath(dest):
                shutil.copy2(built, dest)

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

    python_requires='>=3.9,<3.14',
    scripts=['bin/uls-calc', 'bin/uls-scan', 'bin/uls-nest', 'bin/uls-models', 'bin/uls-scan2D'],

    ext_modules=[quadpack_ext],
    cmdclass={'build_ext': build_ext},
)
