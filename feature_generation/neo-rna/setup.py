# -*- coding: utf-8 -*-

"""
Python Package - Setup
================
"""

from __future__ import print_function

import os
import sys
import imp

from setuptools import setup, find_packages

# ----------------------------------
# region Env. Setup

# Add the current directory to the module search path.
sys.path.insert(0, os.path.abspath('.'))


# Import `metadata`
helper = imp.load_source('project_helper', 'project_helper.py')
metadata = imp.load_source('metadata', os.path.join(helper.CODE_DIRECTORY, 'metadata.py'))

# endregion


# -----------------------------
# `setuptools` setup()

setup_arguments = dict(
    # ----- Basics --------
    name=metadata.package,
    version=metadata.version,

    keywords=metadata.keywords,
    description=metadata.description,
    long_description=helper.read_file('README.rst'),


    author=metadata.authors[0],
    author_email=metadata.emails[0],
    maintainer=metadata.authors[0],
    maintainer_email=metadata.emails[0],

    url=metadata.url,

    # Classifiers
    # :link: https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',

        'License :: OSI Approved :: MIT License',

        'Environment :: Console',

        'Natural Language :: English',
        'Operating System :: OS Independent',

        'Programming Language :: Python :: 2.7',

        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',

        'Topic :: Documentation',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: System :: Installation/Setup',
        'Topic :: System :: Software Distribution',

        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
        'Topic :: Scientific/Engineering :: Medical Science Apps.'
    ],

    packages=find_packages(exclude=(helper.TESTS_DIRECTORY,)),

    install_requires=[
        # Module Dependencies
    ],

    # No build with egg
    zip_safe=False,

    # Entry Point
    entry_points={
        'console_scripts': [
            'neoRNA = neoRNA.main:entry_point'
        ]
    }
)


def main():
    setup(**setup_arguments)


if __name__ == '__main__':
    main()
