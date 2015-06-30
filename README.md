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

6. Comments after `#ifdef`s throw a "extra tokens" compiler warning, so I'm moving the comments to the previous line just to avoid possible problems.

7. Encountering errors in the `read_fort_dg` file; it looks like the "key" value was not property initialized. on line 419 of `read_fort_dg.F90` in the master branch:
`fortdg(1)%key = empty`
should be
`fortdg(i)%key = empty`

But more than this, I think that the string literal "empty" should be used rather than an undefined character sequence.

## List of Data Structures 

SIZES
DG
GLOBAL
WIND
NodalAttributes

## Converting MODULE variables into user defined types

I'm starting with the `SIZES` module, since it's the smallest.  `PARAMETER` variables can't be inside the user defined type, which is fine since they are static. Nothing needs to be done inside the source code that use any of those variables. The rest of the variables are placed inside a block that starts with `type sizes_type`.  No `SAVE`s are allowed.  No `TARGET`s are allowed either, which only affects one variable in this module. The target property seems to be used in reading the keyword style fort.dg file, so I've modified the `read_fort_dg.f95` subroutine, adding a temporary variable `layers_target` and then later assigning that value to `layers`. I think we will need to do this for all of the variables read in via the fort.dg keyword style files.

I've created a list of variables inside `sizes_type` mostly by hand. I've also created a bash script that uses `sed` to prepend `s%` to the beginning of each variable in the list. Each subroutine that uses members of `sizes_type` will need to be modified and passed the instance of `sizes_type` that we will call `s`. 

This process will definitely produce some errors, and need to be manually checked.

This is very useful for interactively replacing a variables in a large number of files: <http://ergoemacs.org/emacs/find_replace_inter.html>

Using this sed command: ` sed 's:\(\bVAR\b\):g%\1:gi'` works fairly well for replacing things VAR with g%VAR. But watch out for repeated "g%"s 

Currently, the procedure I am using involves making a list of variables contained in the user-defined type of the module, for example `src/dg.list`. I have created a script, `scripts/sed_vars.bash` which takes 5 command line arguments, the last two of which are optional: `source_code_filename`, `variable_names_list`, `module_prefix`, `start_line`, and `end_line`. This uses the sed command above to replace all instances of `VAR` in `variable_names_list` with `module_prefix%VAR` in `source_code_filename` between `start_line` and `end_line`.  In instances where a single file contains multiple subroutines, multiple calls to this script can be made with different line numbers specified in order to only replace variables in that subroutine.  In instances where a `USE module_name, ONLY: VAR1, VAR2, ...` is used, a separate list of variables (`VAR1 VAR2 ...`) is made so that only those variables are replaced in that subroutine. Source code files are then manually checked using emacs `ediff` mode.  So far, with module `DG`, the only problems coming from using the sed script is replacing `#ifdef VAR` with `#ifdef dg%VAR`, which needs to be fixed manually.  

The next step is to go to each subroutine using the module and add two things: the definition of the global type `type (global_type) :: global_here`, and adding `global_here` in the subroutine call.

# Module GLOBAL

Will use prefix global_here
global.list

List of variables is `global.list`, script to replace variables is `global_file.list`.

