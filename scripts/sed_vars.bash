#!/bin/bash

# check number of arguments
if [ $# != 2 ]; then
    echo "Usage: $0 source_code_filename variable_list"
    exit 1
fi

sourcecode="$1"
varlist="$2"

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
touch sed_vars_temp_file.txt

tempfile1="sed_vars_temp_file1.txt"
tempfile2="sed_vars_temp_file2.txt"

cp $sourcecode $tempfile1
while read -r var
do
    #echo "$var"
    sed "s/$var/s%$var/g" $tempfile1 > $tempfile2
    cp $tempfile2 $tempfile1
done < "$varlist"

cp $tempfile1 ${sourcecode}.replaced