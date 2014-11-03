import os
import sys
from setuptools import setup

version_py = os.path.join(os.path.dirname(__file__), 'barseqtools', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"','')
long_description = """
``barseqtools`` provides tools for working with barseq data of type indicated in README'
"""

with open("requirements.txt", "r") as f:
    install_requires = [x.strip() for x in f.readlines()]

setup(
        name="barseqtools",
        version=version,
        install_requires=install_requires,
        requires = ['python (>=2.7, <3.0)'],
        packages=['barseqtools',
                  'barseqtools.scripts'],
        author="John Urban (barseqtools)",
        description='Tools for working with barseq data',
        long_description=long_description,
        url="https://github.com/JohnUrban/barseqtools",
        package_dir = {'barseqtools': "barseqtools"},
        package_data = {'barseqtools': []},
        zip_safe = False,
        include_package_data=True,
        #scripts = ['barseqtools/scripts/barseqtools'],
        entry_points = {
            'console_scripts' : [
                 'barseqtools = barseqtools.barseqtools_main:main', 
            ],
        },  
        author_email="mr.john.urban@gmail.com",
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Topic :: Scientific/Engineering :: Bio-Informatics']
    )
