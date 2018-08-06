from setuptools import setup
import os.path

from distutils.command.build_py import build_py

from shapestats import __version__


setup(name='shapestats', # name of package
      version=__version__,
      description='tools & methods to measure shape regularity', #short <80chr description
      url='https://github.com/ljwolf/shapestats', #github repo
      maintainer='Levi John Wolf',
      maintainer_email='levi.john.wolf@gmail.com',
      tests_require=['pytest'],
      keywords='spatial statistics',
      classifiers=[
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: GIS',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6'
        ],
      license='MIT',
      packages=['shapestats'], # add your package name here as a string
      install_requires=['scipy','shapely','libpysal'],
      zip_safe=False,
      cmdclass = {'build.py':build_py})
