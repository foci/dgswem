import subprocess 
import os
import pprint

# This will run all error cases for a convergence study.  
#   - Right now the format of the error.inp file must follow
#     more strict formatting critiera than required by the 
#     error code:
#       1) The file should begin with all case information
#          commented.
#       2) Case information should be written on consecutive 
#          lines with spaces between cases.  No commented lines
#          between case information is allowed
       


nlines = 14         # number of lines per case
file = 'error.inp'  # input file name
exe = 'error'       # executable name



f = open(file,'r+')
content = f.read().splitlines()
f.close()
#pprint.pprint(content)
lines = [x.lstrip().split() for x in content]

# determine number of cases
start_case = 0
ncases = 0
case_ind = []
for n,line in enumerate(lines):
  
  # search for start of case
  if start_case == 0:  
    if not line :  # ignore blank line
      pass
    elif line[0][0]  == '!':  
      if os.path.exists(line[0][1:]): # check if line is a directory           # WARNING: if directory is supposed to exist and 
        ncases = ncases+1             # if it is, it is the start of a case    # doesn't, this wil mess up everything afterward
        start_case = 1                # flag that case has started
        i = 1                         # initalize counter for case lines
        case_ind.append(n)            # keep track of where case begins
        print n
        print line[0]
      else: 
        pass 
    else: 
      print "Format error: all cases should begin commented"
      raise SystemExit(0)
  else:
    if line and line[0][0] == '!':
      i = i+1    
      print line[0]
    else:
      print "Format error: case information needs to be on consecutive lines"
      print "Possible directory exsitence error"
      raise SystemExit(0)

    if i == nlines:  # iterate through case lines until last line of case
      start_case = 0  # flag to search for new case
      print "\n"

   
#print ncases
#print case_ind  







skip = False
# Uncomment and run each case
for n in range(0,ncases):

  ind = case_ind[n]  # look up line where case begins


  print "######################################################"
  print "Case " + str(n+1) + ":"
  for i in range(0,nlines):            # print current case
    print content[ind+i].lstrip()[1:]
  print "######################################################"

  print " "

  print "c) Continue"
  print "s) Skip"
  print "e) Exit"
  option = raw_input("Choose option: ") # pause between each run
  while option != "c" and option != "s" and option != "e":
    option = raw_input("Choose option: ")
  print " "

  if option == "c":
    pass
  elif option == "s":
    skip = True
  elif option == "e":
    break






  if skip == False:

    print "Running case " + str(n+1) + " starting on line " + str(ind)
    print " "

    f = open(file,'r+')

    f.write('\n'.join(content[0:ind])+'\n') # write (commented) lines before case
    #print '\n'.join(content[0:ind])+'\n'
  
    for i in range(0,nlines):                   # write uncommented lines for current case
      f.write(content[ind+i].lstrip()[1:])
      if i < nlines-1:
        f.write('\n')

    f.close()



    try:   # run the error code and display otput
      cmd = ['./'+exe]
      output = subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()[0]
      print output
    except:
      break

  else:

    print "Skipping case " + str(n+1) + " starting on line " + str(ind)
    print " "

    skip = False


  
  


# Write back original file
f = open(file,'r+')
f.write('\n'.join(content))
f.close()

