#!/bin/sh
doxygen dox.config
mv html/search/* search
mv html/*.* .
rm -r html
#python colors.py doxygen.css
#python colors.py tabs.css
#python colors.py search/search.css
