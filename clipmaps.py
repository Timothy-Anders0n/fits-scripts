from astropy.io import fits

maplist = ["HI/lbmap_loc-inter.fits",
           "HI/lbmap_Outer.fits",
           "HI/lbmap_Sgr-Perseus.fits",
           "HI/lbmap_tanSct.fits",
           "CO/comap_loc-inter.fits",
           "CO/comap_Outer.fits",
           "CO/comap_Sgr-Perseus.fits",
           "CO/comap_tanSct.fits",		
           "Dust2/dust_residuals_loc-inter.fits",
           "Dust2/dust_residuals_Outer.fits",
           "Dust2/dust_residuals_Sgr-Perseus.fits",
           "Dust2/dust_residuals_tanSct.fits"]

for filename in maplist:
    hdu = fits.open(filename)[0]
    hdu.data[hdu.data<0.]=0.
    hdu.writeto(filename[:-5]+"_clip.fits")
