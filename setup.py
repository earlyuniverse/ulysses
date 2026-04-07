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
        # build_temp is where temporary files go
        build_temp = pathlib.Path(self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)
        
        # extdir is where the final compiled library goes
        extdir = pathlib.Path(self.get_ext_fullpath(ext.name))
        extdir.mkdir(parents=True, exist_ok=True)
        
        config = 'Debug' if self.debug else 'Release'
        
        # We point cmake to the subdirectory where CMakeLists.txt lives
        cmake_args = [
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + str(extdir.parent.absolute()),
            '-DCMAKE_BUILD_TYPE=' + config
        ]
        
        build_args = ['--config', config]

        os.chdir(str(build_temp))
        # Call cmake pointing to the sourcedir (NumbaQuadpack folder)
        self.spawn(['cmake', ext.sourcedir] + cmake_args)
        if not self.dry_run:
            self.spawn(['cmake', '--build', '.'] + build_args)
        
        # Return to the original directory so setuptools doesn't get lost
        os.chdir(str(pathlib.Path(__file__).parent.absolute()))

# --- Combined Setup ---

setup(
    name='ulysses',
    version='2.0.2',
    description='ULYSSES: Universal LeptogeneSiS Equation Solver',
    url='https://github.com/earlyuniverse/ulysses',
    author='Alessandro Granelli, Christopher Leslie, Yuber Perez Gonzalez, Brian Shuve, Holger Schulz, Jessica Turner, Rosie Walker',
    author_email='jessicaturner.5390@gmail.com',
    
    # find_packages will now pick up both 'ulysses' and 'NumbaQuadpack'
    packages=find_packages(),
    include_package_data=True,
    
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'pymultinest',
        'progressbar',
        'pandas',
        'python-dateutil',
        'termcolor',
        'tqdm',  
        'numba',
        'mpmath',
        'mpltern',
        'multiprocess'
    ],
    
    python_requires='>3.6.0',
    scripts=['bin/uls-calc', 'bin/uls-scan', 'bin/uls-nest', 'bin/uls-models', 'bin/uls-scan2D'],
    
    # This triggers the cmake build for the NumbaQuadpack subfolder
    ext_modules=[CMakeExtension('NumbaQuadpack.cquadpack', sourcedir='NumbaQuadpack')],
    cmdclass={'build_ext': build_ext}
)