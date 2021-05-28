####################################################################################
#Disk
length = 250.                    #AU (size of the cube); length > 2*Rout
Rin = 5.00                       #AU (inner radius of the disk)
Rout = 120.                      #AU (outer radius of the disk)
Rc = 60.                         #AU (transition radius of the disk)
inc = 0.                        #deg (disk inclination between 0 and 90 degrees)
distance = 120.                  #pc  (disk distance)

#Star                            (This is only important if the temperature is computed by the model)
Tstar = 4500.                     #K (Star effective temperature)
Mdot = 1.e-7                     #Solar masses / year
Mstar = 1.                       #Solar masses
Rstar = 6.96e11                  #cm
phi_angle = 0.05                 #rad


#Convolution
beam = 0.1                       #arcsec

####################################################################################
#Radiative transfer
wl = 0.13                        #cm (Observed wavelenght) [it should exist. Options: wl = 0.087, 0.13, 0.3, 0.7]

####################################################################################
#What to do
wr_files = False                 #Write fits files
add_Star = False                 #Add star to the images? Problem with intensity scale
add_back = True                 #Add the 2.7 background intensity?
scattering = True                #Scattering ON or OFF?

plot_opacity = True
plot_sky = True                  #MAP (pdf)
plot_Temperature = True          #MAP of the temperature structure
intensity_log = True             #color map (intensity) in log scale?
opacity_log = True               #color map (opacity) in log scale?

####################################################################################
#Extra parameters (do not change)
nu = clight/wl                    #from wl to nu
inc = inc*np.pi/180.              #from deg to rad
Mdot = Mdot/year*Msun             #from solar masses/year to gram/sec
Mstar = Mstar*Msun                #from solar masses to gram
