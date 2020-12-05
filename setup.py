# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 14:30:27 2018

@author: tadahaya
"""

from setuptools import setup,find_packages
from setuptools.command.test import test as TestCommand
import sys

__version__ = "2.0.0" # major.minor[.patch[.sub]] 

sys.path.append('.\\main')
sys.path.append('.\\tests')

with open('README.rst') as f:
    readme = f.read()

with open('LICENSE.txt') as f:
    license_txt = f.read()

class PyTest(TestCommand):
    user_options = [("pytest-args=", "a", "Arguments to pass to pytest")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = ""

    def run_tests(self):
        import shlex
        # import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(shlex.split(self.pytest_args))
        sys.exit(errno)

setup(
    name='olsapy',
    version=__version__,
    description="a package for orthogonal linear separation analysis (OLSA)",
    long_description=readme,
    url='https://github.com/mizuno-group/OLSApy',
    author='tadahaya mizuno',
    author_email='tadahaya@gmail.com',
    license=license_txt,
    classifiers=['Development Status :: 5 - Production/Stable'
                 ,'Environment :: Win32 (MS Windows)'
                 ,'Framework :: IPython'
                 ,'Intended Audience :: Science/Research'
                 ,'Operating System :: Microsoft :: Windows :: Windows 10'
                 ,'Programming Language :: Python :: 3.6'
                 ,'Topic :: Scientific/Engineering :: Bio-Informatics'],
    keywords=['olsa','usspca','profiling','omics','bioinformatics','profile data','transcriptome'],
    install_requires=["numpy","scipy","pandas"],
    include_package_data=True,
    package_data={'olsapy':['*.txt','*.ignore']},
    packages=find_packages(exclude=('tests', 'docs')),
    tests_require=['pytest'],
    cmdclass={"test": PyTest},
    zip_safe=False
    )