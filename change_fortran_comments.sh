#!/bin/sh

# This tiny script changes any "c" or "C" characters in the first line
# of a fortran code into !s

if [ $# != 1 ]; then
    echo "Usage: $0 filename"
else
    if [ -f $1 ]; then
	echo "Changing comment characters for \"$1\""
	sed 's/^c/!/g' $1 > /tmp/temp
	sed 's/^C/!/g' /tmp/temp > /tmp/temp2
	mv /tmp/temp2 $1
    else
	echo "input file \"$1\" not found"
    fi
fi