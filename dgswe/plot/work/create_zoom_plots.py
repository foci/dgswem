import os
import sys
import shutil
import glob
import subprocess

def write_file(name,content):    
      
      
  if name != '':
    f = open(name,'w')   
    for line in content:
      value = line['value']
      comment = line['comment']    
      spaces = 101 - len(value)    
      f.write(value + spaces*' ' + comment)
    f.close()   
    
    

plot_exe_direc = '/home/sbrus/Codes/dgswe/plot/work/'
plotting_direc = '/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri_x4/adcirc/ESL0/'

zoom_list = [
             ['2.95e5,3.12e5,3.283e6,3.298e6','upper'],  
             ['3.2e5,3.4e5,3.24e6,3.26e6','inlet'],
             ['2.84e5,3.00e5,3.21e6,3.233e6','left inlet']
             ]
             
plot_order = '15'             
             
             

shutil.copy(plot_exe_direc+'plot',plotting_direc)


for zoom in zoom_list :  
  
  zoom_box  = zoom[0]
  zoom_name = zoom[1]
  zoom_name_ = "_".join(zoom_name.split())
  
  save_direc = plotting_direc+'plots/'+zoom_name_+'/'
  if not os.path.exists(save_direc):
    os.makedirs(save_direc)
  
  
  base_cscale_direc = '/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri_x64/p3/ctp3/hbp3/plots/'+zoom_name_+'/'
  shutil.copy(base_cscale_direc+'vel.cscale',plotting_direc)
  shutil.copy(base_cscale_direc+'zeta.cscale',plotting_direc)

  plot_input = [{'value':'./'                                             , 'comment':'!'+'dgswe.inp file path \n'},
                {'value':'/scratch365/sbrus/,/home/sbrus/data-drive/'     , 'comment':'!'+'dgswe.inp path substitution \n'},
                {'value':'1,0,1'                                          , 'comment':'!'+'zeta plot options \n'},
                {'value':'1,0,1'                                          , 'comment':'!'+'velocity plot options \n'},
                {'value':'1,0,1'                                          , 'comment':'!'+'bathymetry plot options \n'},
                {'value':'1'                                              , 'comment':'!'+'mesh plot option \n'},
                {'value':'off'                                            , 'comment':'!'+'mesh plot element labels \n'},
                {'value':'off'                                            , 'comment':'!'+'mesh plot node lables \n'},
                {'value': plot_order                                      , 'comment':'!'+'order of nodal set for plotting straight elements \n'},
                {'value': plot_order                                      , 'comment':'!'+'order of nodal set for plotting curved elements \n'},
                {'value': zoom_box                                        , 'comment':'!'+zoom_name+' zoom box \n'},
                {'value':'7'                                              , 'comment':'!'+'figure width \n'},
                {'value':'1,48'                                           , 'comment':'!'+'timesnap range \n'},
                {'value':'/home/sbrus/Codes/dgswe/plot/work/default2.cmap', 'comment':'!'+'colormap path \n'},
                {'value':'file'                                           , 'comment':'!'+'zeta colorscale limits \n'},      
                {'value':'file'                                           , 'comment':'!'+'velocity colorscale limits \n'},    
                {'value':'11'                                             , 'comment':'!'+'font size \n'},                            
                {'value':'10'                                             , 'comment':'!'+'number of x ticks \n'},              
                {'value':'auto'                                           , 'comment':'!'+'number of y ticks \n'},         
                {'value':'10'                                             , 'comment':'!'+'number of c ticks \n'},      
                {'value':'3'                                              , 'comment':'!'+'number of x decimal places \n'},   
                {'value':'3'                                              , 'comment':'!'+'number of y decimal places \n'},    
                {'value':'3'                                              , 'comment':'!'+'number of c decimal places \n'},  
                {'value':'3'                                              , 'comment':'!'+'number of t decimal places \n'},      
                {'value':'jpg,1'                                          , 'comment':'!'+'additional file format, remove ps file option \n'},  
                {'value':'250'                                            , 'comment':'!'+'raster density \n'},   
                {'value':'0'                                              , 'comment':'!'+'movie flag \n'}        
                ] 
                
  write_file(plotting_direc+'plot.inp' ,plot_input)              

  os.chdir(plotting_direc)
  output = subprocess.Popen('./plot', shell=True, stdout=subprocess.PIPE)
  
  while output.poll() is None:
    out = output.stdout.readline()
    if out:
      print out.strip()
      
    
  subprocess.Popen('mv *.jpg '+save_direc, shell=True, stdout=subprocess.PIPE).communicate()[0]  
  subprocess.Popen('mv *.cscale '+save_direc, shell=True, stdout=subprocess.PIPE).communicate()[0] 
  subprocess.Popen('mv plot.inp '+save_direc, shell=True, stdout=subprocess.PIPE).communicate()[0]   
  subprocess.Popen('cp plot '+save_direc, shell=True, stdout=subprocess.PIPE).communicate()[0]    
  
