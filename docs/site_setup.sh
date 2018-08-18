#!/bin/sh
doxygen dox.config
mv html/search/* search
mv html/*.* .
rm -r html