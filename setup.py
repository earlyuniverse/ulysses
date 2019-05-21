from setuptools import setup, find_packages
setup(
  name = 'ulysses',
  version = '1.0',
  description = 'ULYSSES: Universal LeptogeneSiS Equation Solver',
  url = 'https://github.com/iamholger/ulysses',
  author = 'Kris Moffat, Jessica Turner, Holger Schulz',
  author_email = 'hschulz@fnal.gov',
  packages = find_packages(),
  include_package_data = True,
  install_requires = [
    'numpy',
    'scipy',
    'odeintw',
    'numba',
    'pymultinest'
  ],
  scripts=['bin/uls-calc', 'bin/uls-scan', 'bin/uls-nest'],
  extras_require = {
  },
  entry_points = {
  },
  dependency_links = [
  ]
)
