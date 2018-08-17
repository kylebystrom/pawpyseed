#!/bin/sh
doxygen dox.config
mv html/* .
rm -r html