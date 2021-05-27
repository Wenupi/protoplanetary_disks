from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve

std = beam*distance                       #from arcsec to au
print ('The beam convolution is =', std, 'au')
std = std/2./np.sqrt(2.*np.log(2.))       #from FWHM (beam) to standard deviation
std = std*nx/length                       #beam in pixels

kernel = Gaussian2DKernel(std)

conv_Bright = convolve(Bright,kernel)
conv_TB = hP*nu/kB/np.log(norm/conv_Bright/Jy + 1.)


#The flux should be the same ( ... the sum is approx the same ...)
#conv_Bright = conv_Bright/np.sum(conv_Bright)*np.sum(Bright)
#conv_TB = conv_TB/np.sum(conv_TB)*np.sum(TB)
