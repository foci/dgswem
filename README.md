dgswem release 11.13
=========
Discontinuous Galerkin Shallow Water Equation Model

## Changes made

1. Changed all fortran source code file endings to .f95

2. Changing comment characters (first line = C or c) to !
can do this using `sed 's/^C/!/g' file` and `sed 's/^c/!/g' file`
remove line continuation characters from the start of one line to the end of the previous line.  I created a simple shell script to do this automatically.

3. Replacing line continuation characters from the begining of the continued line to the end of the previous line. Inserting new lines in emacs search and replace is useful (see <http://jeremy.zawodny.com/blog/archives/008872.html> for more information)

In emacs: 
``
search-replace "ctrl-q ctrl-j <five spaces> &" with "& ctrl-q ctrl-j"
search-replace "ctrl-q ctrl-j <five spaces> $" with "& ctrl-q ctrl-j"
``
replaces most of the line-continuation characters. This fails when the previous 
line (the line that needs to be continued) is commented... either the entire line is a comment, or ends in a comment.  

4. "+" characters inside format statements don't seem to be compatible. This seems to just be a line continuation character for FORMAT statements, so using the same search-replace method above but replacing & with + works.

5. Added the following flags to gfortran compiler `-ffixed-line-length-none -ffree-line-length-none`

6. comments after 



## List of Data Structures 
