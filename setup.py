from setuptools import setup
from distutils.util import convert_path

exec(open(convert_path('filters/version.py')).read())

with open('README.md') as f:
    long_description = f.read()

setup(name='filters',
      version=__version__,
      description='Filter response curves',
      long_description=long_description,
      url='https://github.com/adwasser/filters',
      download_url='https://github.com/adwasser/filters/tarball/' + __version__,
      author='Asher Wasserman',
      author_email='adwasser@ucsc.edu',
      license='MIT',
      packages=['filters'],
      package_data={'': ['LICENSE', 'README.md'],
                    'filters': ['version.py', 'data/*']},
      # scripts=['bin/'],
      include_package_data=True,
      install_requires=['numpy', 'scipy', 'matplotlib', 'astropy'],
      zip_safe=False)
