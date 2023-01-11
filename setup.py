from setuptools import setup, find_packages
setup(
  name = 'ulysses',
  version = '2.0.2',
  description = 'ULYSSES: Universal LeptogeneSiS Equation Solver',
  url = 'https://github.com/earlyuniverse/ulysses',
  author = 'Alessandro Granelli, Christopher Leslie, Yuber Perez Gonzalez, Brian Shuve,  Holger Schulz, Jessica Turner, Rosie Walker',
  author_email = 'jessicaturner.5390@gmail.com',
  packages = find_packages(),
  include_package_data = True,
  install_requires = [
    'numpy',
    'scipy',
    'matplotlib',
    'pymultinest',
    'progressbar',
    'python-dateutil',
    'termcolor',
    'tqdm'  
    # 'numba',
    # 'mpi4py'
  ],
  python_requires='>3.6.0',
  scripts=['bin/uls-calc', 'bin/uls-scan', 'bin/uls-nest', 'bin/uls-models', 'bin/uls-scan2D'],
  extras_require = {
  },
  entry_points = {
  },
  dependency_links = [
  ]
)
