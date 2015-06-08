dgswem release 11.13
=========
Discontinuous Galerkin Shallow Water Equation Model

## Changes made

Changed global.F and dg.F to global.f95 and dg.f95 respectively,
changing comment characters (first line = C or c) to !
can do this using `sed 's/^C/!/g' file` and `sed 's/^c/!/g' file`
remove line continuation characters from the start of one line to the end of the previous line.  emacs macros are useful for this, as is searching for "     &".

Insert a newline into an emacs search and replace with this:
http://jeremy.zawodny.com/blog/archives/008872.html

## List of Data Structures 
