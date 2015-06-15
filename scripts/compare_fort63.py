#!/usr/bin/env python

hpx_fort63_file = "../../dgswem_example/fort.63"
reference_fort63_file = "../../dgswem_example_reference/fort.63"

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
    print diff
    if (diff > max_relative_diff):
        max_relative_diff = diff

print "max_relative_diff = ", max_relative_diff
