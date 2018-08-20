#!/bin/sh
doxygen dox.config
mv html documentation
cat ../README.md > index.md
