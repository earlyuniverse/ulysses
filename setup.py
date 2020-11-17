from setuptools import setup, find_packages
setup(
  name = 'ulysses',
  version = '1.0.1',
  description = 'ULYSSES: Universal LeptogeneSiS Equation Solver',
  url = 'https://github.com/earlyuniverse/ulysses',
  author = 'Kris Moffat, Yuber Perez Gonzalez, Alessandro Granelli, Holger Schulz, Jessica Turner',
  author_email = 'jessicaturner.5390@gmail.com',
  packages = find_packages(),
  include_package_data = True,
  install_requires = [
    'numpy',
    'scipy',
    'matplotlib',
    'pymultinest',
    'progressbar'
    # 'numba',
    # 'mpi4py'
  ],
  python_requires='>3.6.0',
  scripts=['bin/uls-calc', 'bin/uls-scan', 'bin/uls-nest', 'bin/uls-models'],
  extras_require = {
  },
  entry_points = {
  },
  dependency_links = [
  ]
)
