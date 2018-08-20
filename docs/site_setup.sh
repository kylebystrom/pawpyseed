#!/bin/sh
doxygen dox.config
rm -r documentation/*
mv html/* documentation
rm -r html
cat ../README.md > index.md
