#!/usr/bin/env python


exec(open('libraries.py').read())
exec(open('constants.py').read())

####################################################################################################

start_time = time.time()

####################################################################################################
#Input parameters
exec(open('parameters.py').read())
exec(open('equations.py').read())
exec(open('functions.py').read())
print ('**********************************************************')
print ('The disk inclination is', inc*180./np.pi, 'degrees')
print ('The distance to the disk is', distance, 'pc')
print ('The wavelenght is', wl, 'cm, frequency = ', nu/1.e9, 'GHz')

####################################################################################################
#Sky plane
nx = 2**7+1 #+1 para tener un centro
ny = 2**7+1
x_array = np.zeros((nx))
y_array = np.zeros((ny))

for i in range(nx): #Llena la matriz con valores fisicos
    x_array[i] = length*(i/(nx-1.)-0.5)
for j in range(ny):
    y_array[j] = length*(j/(ny-1.)-0.5)

####################################################################################################
#Creates the 3D Source function
print ('- Creating 3D source function ...')
funcion_S, funcion_T = FUN_creates_source_function(x_array,y_array)
####################################################################################################
#2D grid (sky plane)
Bright = np.zeros((ny,nx))
op_depth = np.zeros((ny,nx))
norm = 2.*hP*nu**3/clight**2
print ('- Creating image ...')
for i in range(nx):
#for i in range(50,51):
    print (np.round(i*100./nx,2), '%')
    for j in range(ny):
        z_crit = FUN_limits_mult(x_array[i],y_array[j])
        if ( len(z_crit) == 0): #i.e. No density in the line of sight
            if (add_back == True):
                Bright[j][i] = FUN_BB(nu,2.7)*norm
            else:
                Bright[j][i] = np.nan
            op_depth[j][i] = np.nan
        else:
            Nhalf = int(len(z_crit)/2)
            
            #Creates arrays of Z at each range
            z_array = np.zeros((Nhalf,100)) 
            for iz in range(Nhalf):
                z_array[iz] = np.linspace(z_crit[2*iz],z_crit[2*iz+1],100)
                
            #Saves the vector optdepth at each range
            sol = np.zeros((Nhalf,100)) 
            for iz in range(Nhalf):
                sol[iz] = odeint(FUN_tau,0.,z_array[iz],
                                 args=(x_array[i],y_array[j])).T[0]
                
            #Saves the total optdepth at each range
            opt = np.zeros((Nhalf))
            for iz in range(Nhalf):
                opt[iz] = sol[iz][len(z_array[iz])-1]
                
            #Saves the total optdepth in the line of sight
            op_depth[j][i] =  np.sum(opt)

            #Saves the optdepth functions (to be interpolated)
            optde = []                      
            for iz in range(Nhalf):
                optde.append(interpolate.interp1d(z_array[iz],sol[iz],kind='linear',fill_value='extrapolate')) #quitar extrapolate?
                
            #Saves the total intensity at each range
            sol_I = np.zeros((Nhalf,100)) 
            for iz in range(Nhalf):
                sol_I[iz] = odeint(FUN_intensity,0.,z_array[iz],
                                   args=(x_array[i],y_array[j],optde[iz])).T[0]*norm
                                
            #Sum all the intensities in the line of sight
            #It takes into account all the right-side optical depth regions
            for iz in range(Nhalf):
                f_ext = 1.
                for indice in range(0,iz):
                    f_ext = f_ext*np.exp(-opt[indice])
                Bright[j][i] = Bright[j][i] + sol_I[iz][len(z_array[iz])-1]*f_ext
                
            #Add the background radiation
            if (add_back == True):
                Bright[j][i] = Bright[j][i] + FUN_BB(nu,2.7)*np.exp(-op_depth[j][i])*norm

        ######################################################################################################
        #Agrega la estrella
        ######################################################################################################
        if((i == (nx-1)/2) and (j == (nx-1)/2) and (add_Star == True)):
            if(len(z_crit) == 0):
                opt_star = 0.
            else:
                if(len(z_crit) == 2): #No hay casos si se puede concatenar todo de alguna manera inteligente
                    opt_star = sol[0]
                    z_star = z_array[0]
                else: #En realidad ANTES de llegar a la estrella solo puede ser 2 o 4 (en un disco)
                    sol[1] = sol[1] + sol[0][len(sol[0])-1] #accumulated optical depth
                    opt_star = np.concatenate((sol[0],sol[1]))
                    z_star = np.concatenate((z_array[0],z_array[1]))
                z_star = z_star[np.where(z_star>=0.)] #The star is at z=0
                opt_star = opt_star[np.where(z_star>=0.)]
                opt_star = opt_star[len(opt_star)-1]
        ######################################################################################################
            B_star = FUN_BB(nu,Tstar)*np.exp(-opt_star)*norm    
            Bright[j][i] = Bright[j][i] + B_star
            print ('   - Star added to the image:')
            print ('     tau on the star = ', opt_star)
                
Area  = (x_array[1] - x_array[0])*(y_array[1] - y_array[0])*AU**2
Bright = Bright*Area/(distance*pc)**2/Jy

if (add_Star == True):
    B_star = B_star*Area/(distance*pc)**2/Jy
    print ('- The disk + star flux is = '+str(np.nansum(Bright)*1.e3)+' mJy')
    print ('- The star flux is = '+str(B_star*1.e3)+'mJy')
    print ('- The disk star flux is = '+str((np.nansum(Bright)-B_star)*1.e3)+' mJy')
else:
    print ('- The total disk flux is = '+str(np.nansum(Bright)*1.e3)+' mJy')



####################################################################################################
#From Jy/pix to Kelvin/pix
TB = hP*nu/kB/np.log(norm/Bright/Jy + 1.)

####################################################################################################
#Map of the temperature structure
z_temp = np.linspace(0.,EQ_Height(np.max(x_array)),100)
map_T = np.zeros((len(z_temp),len(x_array)))

for k in range(len(z_temp)):
    for i in range(len(x_array)):
        if (abs(x_array[i]) < Rin or abs(x_array[i]) > Rout or  z_temp[k] > 2.*EQ_Height(abs(x_array[i]))):
            map_T[k,i] = np.nan
        else:
            map_T[k,i] = EQ_temperature(x_array[i],0.,z_temp[k])


####################################################################################################
#Convolves
exec(open('convolution.py').read())
            
####################################################################################################
print ('- Files written:')
if plot_sky ==True:
#if( (plot_opacity == True) or (plot_sky ==True)):
    exec(open('plots.py').read())
if(wr_files ==True):
    exec(open('write_files.py').read())
####################################################################################################
final_time = time.time()
print ('- Total execution time =', (final_time - start_time)/60., 'minutes')
print ('**********************************************************')
####################################################################################################
