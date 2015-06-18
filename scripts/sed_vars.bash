#!/bin/bash

# check number of arguments
#if [ $# != 2 ]; then
#    echo "Usage: $0 source_code_filename variable_list prefix startline endline"
#    exit 1
#fi

sourcecode="$1"
varlist="$2"
prefix="$3"

if [ $# == 3 ]; then
    startline=1
    endline=`wc -l $sourcecode | awk '{print $1}'`
else
    startline="$4"
    endline="$5"
fi


echo "source code file is $sourcecode and variable list is $varlist"

if [ ! -f $2 ]; then
    echo "variable list \"$2\" not found"
    exit 1
fi

if [ ! -f $1 ]; then
    echo "source code file \"$1\" not found"
    exit 1
fi

# loop through cariable names in 
echo "Changing variables for \"$1\""

# create temporary file

tempfile1="sed_vars_temp_file1.txt"
tempfile2="sed_vars_temp_file2.txt"

echo "sourcecode is \"$sourcecode\", varlist is \"$varlist\", startline is \"$startline\", endline is \"$endline\", prefix is \"$prefix\""

cp $sourcecode $tempfile1
while read -r var
do
#    echo "${startline},${endline}s:\(\b${var}\b\):${prefix}%\1:gi"
    sed "${startline},${endline}s:\(\b${var}\b\):${prefix}%\1:gi" $tempfile1 > $tempfile2	
    cp $tempfile2 $tempfile1
done < "$varlist"

cp $tempfile1 ${sourcecode}.replaced