from astropy.io import fits

fits.writeto('Results/fits/Intensity_wl'+str(np.round(wl,2))+'_inc'+str(np.round(inc*180./np.pi,2))+'deg.fits',
             Bright, overwrite=True)
print ('   - Results/fits/Intensity_wl'+str(np.round(wl,2))+'_inc'+str(np.round(inc*180./np.pi,2))+'deg.fits saved!')

fits.writeto('Results/fits/OpticalD_wl'+str(np.round(wl,2))+'_inc'+str(np.round(inc*180./np.pi,2))+'deg.fits',
             op_depth, overwrite=True)
print ('   - Results/fits/OpticalD_wl'+str(np.round(wl,2))+'_inc'+str(np.round(inc*180./np.pi,2))+'deg.fits saved!')
