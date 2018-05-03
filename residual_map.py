from astropy.io import fits
import numpy as np
import aplpy

counts_list = ["bin_300-605/countMap_flat_PSF3.fits",
		"bin_300-605/countMap_flat_PSF02.fits",
		"bin_605-1220/countMap_flat_PSF3.fits",
		"bin_605-1220/countMap_flat_PSF02.fits",
		"bin_1220-2460/countMap_flat_PSF3.fits",
		"bin_1220-2460/countMap_flat_PSF02.fits",
		"bin_2460-4959/countMap_flat_PSF3.fits",
		"bin_2460-4959/countMap_flat_PSF02.fits",
		"bin_4959-10000/countMap_flat_PSF3.fits",
		"bin_4959-10000/countMap_flat_PSF02.fits"]
		
models_list = ["bin_300-605/modelMap_PSF3.fits",
		"bin_300-605/modelMap_PSF02.fits",
		"bin_605-1220/modelMap_PSF3.fits",
		"bin_605-1220/modelMap_PSF02.fits",
		"bin_1220-2460/modelMap_PSF3.fits",
		"bin_1220-2460/modelMap_PSF02.fits",
		"bin_2460-4959/modelMap_PSF3.fits",
		"bin_2460-4959/modelMap_PSF02.fits",
		"bin_4959-10000/modelMap_PSF3.fits",
		"bin_4959-10000/modelMap_PSF02.fits"]

for s in range(len(counts_list)):
    if s==0:
        hdu_counts = fits.open(counts_list[s])[0]
        model = fits.getdata(models_list[s],0)
    else:
        hdu_counts.data += fits.getdata(counts_list[s],0)
        model += fits.getdata(models_list[s],0)
		
hdu_counts.data = (hdu_counts.data - model)/np.sqrt(hdu_counts.data)

hdu_counts.writeto('residuals.fits')

print('\n Currently processing residuals')
fig = aplpy.FITSFigure(hdu_counts)

fig.show_colorscale(cmap='bwr', vmin=-4, vmax=4)

fig.show_regions('selectedreg.reg')

#fig.refresh() # Forces refresh

#fig.set_auto_refresh(True) # Automatic refresh after each command

fig.add_colorbar()
#fig.colorbar.set_location('bottom') # Location
fig.colorbar.set_axis_label_text('Significance ($\sigma$)') # Label

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

fig.save("residuals.png")

print('Script successful')
