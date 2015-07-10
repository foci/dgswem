#!/usr/bin/env python

# ********* PE0000 ***************

hpx_fort63_file = "../../dgswem_example_4_hpx/PE0000/fort.63"
reference_fort63_file = "../../dgswem_example_4_cpp_series/PE0000/fort.63"

with open(hpx_fort63_file) as f:
    hpx = f.readlines()

with open(reference_fort63_file) as f2:
    ref = f2.readlines()

ref_data = []
hpx_data = []

for line in hpx:
    try:
        wse = float( line.rstrip('\n').split()[1] )
    except ValueError:
        continue
    hpx_data.append(wse)
    
for line in ref:
    try:
        wse = float( line.rstrip('\n').split()[1] )
    except ValueError:
        continue
    ref_data.append(wse)

print "len(ref_data) = ", len(ref_data)
print "len(hpx_data) = ", len(hpx_data)

max_relative_diff = 0.0
for x in xrange(len(ref_data)):
    delta = ref_data[x] - hpx_data[x]
    diff = abs(delta/max(ref_data[x],hpx_data[x]))
#    print diff
    if (diff > max_relative_diff):
        max_relative_diff = diff

print "PE0000 max_relative_diff = ", max_relative_diff

#************************************

# ********* PE0001 ***************

hpx_fort63_file = "../../dgswem_example_4_hpx/PE0001/fort.63"
reference_fort63_file = "../../dgswem_example_4_cpp_series/PE0001/fort.63"

with open(hpx_fort63_file) as f:
    hpx = f.readlines()

with open(reference_fort63_file) as f2:
    ref = f2.readlines()

ref_data = []
hpx_data = []

for line in hpx:
    try:
        wse = float( line.rstrip('\n').split()[1] )
    except ValueError:
        continue
    hpx_data.append(wse)
    
for line in ref:
    try:
        wse = float( line.rstrip('\n').split()[1] )
    except ValueError:
        continue
    ref_data.append(wse)

print "len(ref_data) = ", len(ref_data)
print "len(hpx_data) = ", len(hpx_data)

max_relative_diff = 0.0
for x in xrange(len(ref_data)):
    delta = ref_data[x] - hpx_data[x]
    diff = abs(delta/max(ref_data[x],hpx_data[x]))
#    print diff
    if (diff > max_relative_diff):
        max_relative_diff = diff

print "PE0001 max_relative_diff = ", max_relative_diff

#************************************

# ********* PE0002 ***************

hpx_fort63_file = "../../dgswem_example_4_hpx/PE0002/fort.63"
reference_fort63_file = "../../dgswem_example_4_cpp_series/PE0002/fort.63"

with open(hpx_fort63_file) as f:
    hpx = f.readlines()

with open(reference_fort63_file) as f2:
    ref = f2.readlines()

ref_data = []
hpx_data = []

for line in hpx:
    try:
        wse = float( line.rstrip('\n').split()[1] )
    except ValueError:
        continue
    hpx_data.append(wse)
    
for line in ref:
    try:
        wse = float( line.rstrip('\n').split()[1] )
    except ValueError:
        continue
    ref_data.append(wse)

print "len(ref_data) = ", len(ref_data)
print "len(hpx_data) = ", len(hpx_data)

max_relative_diff = 0.0
for x in xrange(len(ref_data)):
    delta = ref_data[x] - hpx_data[x]
    if ((ref_data[x] == 0.0) | (hpx_data[x] == 0.0)):
        diff = 0.0
    else:
        diff = abs(delta/max(ref_data[x],hpx_data[x]))
#    print "diff_here = ", diff
    if (diff > max_relative_diff):
        max_relative_diff = diff

print "PE0002 max_relative_diff = ", max_relative_diff

#************************************

# ********* PE0003 ***************

hpx_fort63_file = "../../dgswem_example_4_hpx/PE0003/fort.63"
reference_fort63_file = "../../dgswem_example_4_cpp_series/PE0003/fort.63"

with open(hpx_fort63_file) as f:
    hpx = f.readlines()

with open(reference_fort63_file) as f2:
    ref = f2.readlines()

ref_data = []
hpx_data = []

for line in hpx:
    try:
        wse = float( line.rstrip('\n').split()[1] )
    except ValueError:
        continue
    hpx_data.append(wse)
    
for line in ref:
    try:
        wse = float( line.rstrip('\n').split()[1] )
    except ValueError:
        continue
    ref_data.append(wse)

print "len(ref_data) = ", len(ref_data)
print "len(hpx_data) = ", len(hpx_data)

max_relative_diff = 0.0
for x in xrange(len(ref_data)):
    delta = ref_data[x] - hpx_data[x]
    if ((ref_data[x] == 0.0) | (hpx_data[x] == 0.0)):
        diff = 0.0
    else:
        diff = abs(delta/max(ref_data[x],hpx_data[x]))
#    print diff
    if (diff > max_relative_diff):
        max_relative_diff = diff

print "PE0003 max_relative_diff = ", max_relative_diff

#************************************
