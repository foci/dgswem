import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib import interactive
import numpy as np
import os
import curses

interactive(True)

frmt = 'jpg'
sol = 'vel'

direc_ls = [
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_tri_x64/p3/ctp3/hbp3/plots/upper/',
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_tri/p3/ctp3/hbp3/plots/upper/',
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_tri/p3/ctp1/hbp1/plots/upper/',
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_quad/p3/ctp3/hbp3/plots/upper/',  
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_quad/p3/ctp1/hbp1/plots/upper/',        
            
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_tri_x64/p3/ctp3/hbp3/plots/upper/',
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_tri/p3/ctp3/hbp3/plots/upper/',
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_tri/p2/ctp2/hbp2/plots/upper/',
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_tri/p1/ctp2/hbp1/plots/upper/',  
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_tri_x4/p3/ctp3/hbp3/plots/upper/',
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_tri_x4/p2/ctp2/hbp2/plots/upper/',
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_tri_x4/p1/ctp2/hbp1/plots/upper/',            
            
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_tri_x64/p3/ctp3/hbp3/plots/inlet/',
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_tri_x4/p3/ctp3/hbp3/plots/inlet/',
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_tri/p3/ctp3/hbp3/plots/inlet/',                      
            
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_tri_x64/p3/ctp3/hbp3/plots/inlet/',
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_tri/p3/ctp3/hbp3/plots/inlet/',
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_tri/p3/ctp1/hbp1/plots/inlet/',            
            
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_tri_x64/p3/ctp3/hbp3/plots/left_inlet/',
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_tri/p3/ctp3/hbp3/plots/left_inlet/',
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_tri/p3/ctp1/hbp1/plots/left_inlet/'      
            
            #'/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri_x64/p3/ctp3/hbp3/plots/upper/',
            #'/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri_x16/p3/ctp3/hbp3/plots/upper/',
            #'/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri_x4/p3/ctp3/hbp3/plots/upper/',
            #'/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri/p3/ctp3/hbp3/plots/upper/',  
            #'/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_quad_x64/p3/ctp3/hbp3/plots/upper/',
            #'/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_quad_x16/p3/ctp3/hbp3/plots/upper/',
            #'/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_quad_x4/p3/ctp3/hbp3/plots/upper/',
            #'/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_quad/p3/ctp3/hbp3/plots/upper/',          
            
            #'/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri_x64/adcirc/ESL0/plots/upper/',
            #'/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri_x64/p3/ctp3/hbp3/plots/upper/',
            #'/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_quad_x64/p3/ctp3/hbp3/plots/upper/',       
            
            #'/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri_x64/adcirc/ESL0/plots/inlet/',
            #'/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri_x64/p3/ctp3/hbp3/plots/inlet/',
            #'/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_quad_x64/p3/ctp3/hbp3/plots/inlet/',    
            
            #'/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_quad_x64/p3/ctp3/hbp3/plots/inlet/',  
            #'/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri/p3/ctp3/hbp3/plots/inlet/',   
            #'/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri/p2/ctp2/hbp2/plots/inlet/',              
            #'/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri/p1/ctp2/hbp1/plots/inlet/',             
            #'/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri/adcirc/ESL0/plots/inlet/',        
            
            '/home/sbrus/data-drive/galveston_SL18_tides/galveston_SL18/p1/ctp2/hbp1/',
            '/home/sbrus/data-drive/galveston_SL18_tides/galveston_SL18/p1/ctp1/hbp1/',
            ]         
            
#offset = [0,0,0]         
#offset = [1,1,1,1,0] #account for adcirc not writing the initial condition

offset = []
for direc in direc_ls:
  if 'adcirc' in direc:
    offset.append(0)
  else:
    offset.append(1)
    
    
keep_going = True  
snap = 1
while keep_going:          
  
  for i,direc in enumerate(direc_ls):         
  
    snap_str = "%04d" % (snap+offset[i])  

    im_file = direc + sol + '_' + snap_str + '.' + frmt 
    if not os.path.isfile(im_file):
      keep_going = False
      break
  
    fig = plt.figure(i)
    fig.canvas.set_window_title(direc)
    img=mpimg.imread(im_file)      
    plt.imshow(img)
    plt.tight_layout()
    plt.axis('off')


  inp = raw_input("Next image: >, Previous image: <  ")   
  while inp != ',' and inp != '.':
    inp = raw_input("Next image: >, Previous image: <  ") 
    
  if inp == '.':
    snap = snap + 1
  elif inp == ',' and snap > 1:
    snap = snap - 1

