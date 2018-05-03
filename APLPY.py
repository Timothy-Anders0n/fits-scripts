import matplotlib
matplotlib.use('Agg')
import aplpy
import glob
import os
from astropy import units as u
from multiprocessing import Pool

processes = 11

for filepath in glob.iglob('./**/*.fits', recursive=True):
    try:
        print('\n Currently processing: ' + filepath + '\n')
        fig = aplpy.FITSFigure(filepath)
        if '_clip.fits' in filename:
            fig.show_colorscale(cmap='inferno', vmin=0)
        else:
            fig.show_colorscale(cmap='inferno')

        #fig.show_regions('myregions.reg')

        #fig.refresh() # Forces refresh

        #fig.set_auto_refresh(True) # Automatic refresh after each command

        fig.add_colorbar()
        #fig.colorbar.set_location('bottom') # Location
        if './HI/' in filepath:    
            fig.colorbar.set_axis_label_text('HI Column Density ($10^{20} \; HI.cm^{-2}$)') # Label
        elif './CO/' in filepath:
            fig.colorbar.set_axis_label_text('CO Line Intensity ($K.km.s^{-1}$)') # Label
        elif './Dust2/' in filepath:
            fig.colorbar.set_axis_label_text('Dust Opacity (unit of magnitude $mag$)') # Label
        elif 'Count' in filepath:
            fig.colorbar.set_axis_label_text('Photon counts (unit of photon)') # Label
        elif 'res' or 'sig' or 'sqmodel' or 'frac' in filepath:
            fig.colorbar.set_axis_label_text('Significance ($\sigma$)') # Label
        else:
            fig.colorbar.set_axis_label_text('Unregistered Quantity') # Label    

        #fig.add_scalebar()
        #fig.scalebar.set_corner('top right') # Location
        #fig.scalebar.set_frame(True) # Frame
        #fig.scalebar.set_alpha(0.7) # Transparency
        #fig.scalebar.show()
    
        fig.set_xaxis_coord_type('longitude')
        fig.set_yaxis_coord_type('latitude')
    
        fig.axis_labels.show()
        fig.axis_labels.set_xtext('Longitude (GAL)') # Label X axis
        fig.axis_labels.set_ytext('Latitude (GAL)') # Label Y axis
        fig.tick_labels.set_xformat('dd')
        fig.tick_labels.set_yformat('dd')
        fig.tick_labels.set_font(size='large')
        fig.axis_labels.set_font(size='large')
        fig.colorbar.set_font(size='large')


        #fig.set_theme('publication') # Change look of the plot

        fig.save(filepath.replace('.fits','.png'))
        print('end \n \n')

    except:
        pass

print('Script successful')
