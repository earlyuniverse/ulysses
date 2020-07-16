from setuptools import setup, find_packages
setup(
  name = 'ulysses',
  version = '1.0.0',
  description = 'ULYSSES: Universal LeptogeneSiS Equation Solver',
  url = 'https://github.com/earlyuniverse/ulysses',
  author = 'Kris Moffat, Yuber Perez Gonzalez, Alessandro Granelli, Holger Schulz, Jessica Turner',
  author_email = 'jessicaturner.5390@gmail.com',
  packages = find_packages(),
  include_package_data = True,
  install_requires = [
    'numpy',
    'scipy',
    'numba',
    'matplotlib',
    'pymultinest',
    'mpi4py'
  ],
  python_requires='>3.5.2',
  scripts=['bin/uls-calc', 'bin/uls-scan', 'bin/uls-nest', 'bin/uls-models'],
  extras_require = {
  },
  entry_points = {
  },
  dependency_links = [
  ]
)
