# -*- coding: utf-8 -*-
"""Project metadata

Information to describe this project.

NOTE:
    Be sure to "not" include any dependency here, otherwise they will need to be added to
    `setup.py` as well.

    To import it to `setup.py`, use `imp.load_source()` method to directly load it as a module.

"""

# The package name, which is also the "UNIX name" for the project.
package = 'neoRNA'

project = "neoRNA"
project_no_spaces = project.replace(' ', '_')

version = '0.6.6'

keywords = 'NGS pipeline RNA'
description = 'A Python Toolkit for RNA NGS and Analysis.'

authors = ['Lei SHI']
authors_string = ', '.join(authors)
emails = ['foxshee@gmail.com']

license = 'MIT'
copyright = '2017 ' + authors_string

url = 'http://example.com/'
