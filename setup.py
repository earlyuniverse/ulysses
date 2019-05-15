from setuptools import setup, find_packages
setup(
  name = 'leptomts',
  version = '0.3.0',
  description = 'leptMTS',
  url = 'https://github.com/iamholger/leptomts',
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
  scripts=['bin/lepto-get', 'bin/lepto-scan'],
  extras_require = {
  },
  entry_points = {
  },
  dependency_links = [
  ]
)
