import os
import pathlib
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext as build_ext_orig

# --- NumbaQuadpack Custom Build Logic ---

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        super().__init__(name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class build_ext(build_ext_orig):
    def run(self):
        for ext in self.extensions:
            if isinstance(ext, CMakeExtension):
                try:
                    self.build_cmake(ext)
                except Exception as e:
                    print(
                        "\n*** Warning: NumbaQuadpack C extension could not be built ***\n"
                        f"  Reason: {e}\n"
                        "  Models requiring NumbaQuadpack (e.g. Case2,Case3,Case4,CaseS2) will not be available.\n"
                        "  To build manually: cd NumbaQuadpack && cmake . && cmake --build .\n"
                    )
        super().run()

    def build_cmake(self, ext):
        build_temp = pathlib.Path(self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)

        # Build the library directly into the NumbaQuadpack source directory
        # so that (a) editable installs work immediately, and (b) pip copies
        # it to site-packages as package_data on a regular install.
        pkg_dir = pathlib.Path(__file__).parent / "NumbaQuadpack"

        config = 'Debug' if self.debug else 'Release'
        cmake_args = [
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + str(pkg_dir.absolute()),
            '-DCMAKE_BUILD_TYPE=' + config,
        ]
        build_args = ['--config', config]

        os.chdir(str(build_temp))
        self.spawn(['cmake', ext.sourcedir] + cmake_args)
        if not self.dry_run:
            self.spawn(['cmake', '--build', '.'] + build_args)
        os.chdir(str(pathlib.Path(__file__).parent.absolute()))

# --- Combined Setup ---

setup(
    name='ulysses',
    version='2.0.2',
    description='ULYSSES: Universal LeptogeneSiS Equation Solver',
    url='https://github.com/earlyuniverse/ulysses',
    author='Alessandro Granelli, Christopher Leslie, Yuber Perez Gonzalez, Brian Shuve, Holger Schulz, Jessica Turner, Rosie Walker',
    author_email='jessicaturner.5390@gmail.com',
    
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'NumbaQuadpack': ['*.so', '*.dylib', '*.dll'],
    },
    
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
    
    # This triggers the cmake build for the NumbaQuadpack subfolder
    ext_modules=[CMakeExtension('NumbaQuadpack.cquadpack', sourcedir='NumbaQuadpack')],
    cmdclass={'build_ext': build_ext}
)