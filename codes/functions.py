#!/usr/bin/env python

#--------------------------------------------------------------------------------
#Changes the sky coordinates (x,y,z) to the disk coordinates (x_d,y_d,z_d)
#The x axis is the rotation axis
def FUN_rotation(x,y,z):
    x_d = x    
    y_d = y*np.cos(inc) - z*np.sin(inc)
    z_d = y*np.sin(inc) + z*np.cos(inc)
    return x_d,y_d,z_d

#--------------------------------------------------------------------------------
#Radiative transfer equation
def FUN_intensity(I,z,x,y,optde):
    x_d,y_d,z_d = FUN_rotation(x,y,z)
    density = EQ_density(x_d,y_d,z_d)
    amax = EQ_amax(x_d,y_d,z_d)
    opa = function_ext(amax)
    S = funcion_S([z_d,y_d,x_d])
#    print ('x,y,z', x,y,z)
#    print (S, x_d, y_d, z_d)
#    print (optde(z))
    dIdz = -S*opa*density*np.exp(-optde(z)) #z es la variable de integracion (debe ser evaluada en cualquier punto)
    return dIdz

#--------------------------------------------------------------------------------
#Optical depth
def FUN_tau(tt,z,x,y):
    x_d,y_d,z_d = FUN_rotation(x,y,z)
    density = EQ_density(x_d,y_d,z_d)
    amax = EQ_amax(x_d,y_d,z_d)
    opa = function_ext(amax)
    dtau = -opa*density
    return dtau

#--------------------------------------------------------------------------------
def FUN_tau_zaxis(tt,z,x,y):
    x_d,y_d,z_d = x,y,z
    density = EQ_density(x_d,y_d,z_d)
    amax = EQ_amax(x_d,y_d,z_d)
    opa = function_ext(amax)
    dtau = -opa*density
    return dtau

#--------------------------------------------------------------------------------
#Black body radiation
def FUN_BB(nu,T):
#    B = 2.*hP*nu**3/clight**2/( np.exp(hP*nu/kB/T) - 1.)
    B = 1./( np.exp(hP*nu/kB/T) - 1.)
    return B

#--------------------------------------------------------------------------------
def FUN_limits_mult(xx,yy):
    Hout = EQ_Height(Rout)
    lim_z = Rout*np.sin(inc) + 2.*Hout*np.cos(inc) #Based on the geometry of the disk
    lim_y = Rout*np.cos(inc) + 2.*Hout*np.sin(inc)  #Based on the geometry of the disk

    z_arr = np.linspace(1.1*lim_z, -1.1*lim_z, 200)

    z_crit = []

    if ((np.abs(xx) <=Rout) and (np.abs(yy) <= lim_y)):
        xd,yd,zd = FUN_rotation(xx,yy,z_arr)
        crit = np.zeros((len(z_arr)))
        ###############################################################################
        #Funciona pero podria ser optimizado
        ###############################################################################        
        for ii in range(len(z_arr)): #Crea un vector de densidad en la linea de vision
            if (EQ_density(xd,yd[ii],zd[ii]) == 0.):
                crit[ii] = 0
            else:
                crit[ii] = 1
    
        for ii in range(len(z_arr)): #Ve los indices donde cambia de 0 a algun valor, o de algun valor a 0 (fronteras)
            if ( (ii != 0) and (crit[ii] - crit[ii-1] != 0 )):
                z_crit.append(z_arr[ii])
            elif(ii == 0 and crit[0] == 1):
                z_crit.append(z_arr[0])
        ###############################################################################                
                
    return z_crit

#--------------------------------------------------------------------------------
def FUN_creates_source_function(x_array,y_array):
    #Arrays and limits
    Hout = EQ_Height(Rout)
    z_array = np.linspace(-2.*Hout, 2.*Hout, 200)
    Sfunction = np.zeros((len(z_array),len(y_array),len(x_array)))
    Temfunction = np.zeros((len(z_array),len(y_array),len(x_array)))    
    op_depth_p = np.zeros((len(y_array),len(x_array)))

    #Computes the optical depth (perpendicular to the disk midplane)
    for j in range(len(y_array)):
        for i in range(len(x_array)):
            if(x_array[i] == 0. and y_array[j] == 0.):
                Sfunction[:,j,i] = 0.
                Temfunction[:,j,i] = 0.
            else:
                rad = np.sqrt(x_array[i]**2 + y_array[j]**2)
                Hscale = EQ_Height(rad)

                z_integ = np.linspace(2.*Hscale,-2.*Hscale,200)           
                sol = odeint(FUN_tau_zaxis,0.,z_integ,args=(x_array[i],y_array[j])).T[0]
                op_depth_p[j][i] = sol[len(z_integ)-1]

                inter_opt = interpolate.interp1d(z_integ,sol,kind='linear', bounds_error=False,fill_value=0.)
                for k in range(len(z_array)):
                    amax = EQ_amax(x_array[i],y_array[j],z_array[k])
                    albedo = function_alb(amax)
                    ##########Temperature##########
                    Omega2 = Ggrav*Mstar/(rad*AU)**3
                    Teff4 = 3.*Mdot*Omega2/8./np.pi/sigmaB
                    Tacc4 = 3./4.*(7.*inter_opt(abs(z_array[k])) + 2./3.)*Teff4
                    Tirr4 = Tstar**4./4.*(Rstar/rad/AU)**2*np.exp(-7.*inter_opt(abs(z_array[k]))/phi_angle)
                    Temfunction[k,j,i] = (Tacc4 + Tirr4)**(0.25)
                    #Temfunction[k,j,i] = EQ_temperature(x_array[i],y_array[j],z_array[k])
                    ###############################
                    Sfunction[k,j,i] = FUN_BB(nu,Temfunction[k,j,i])*(1.+ albedo*FUN_f(inter_opt(z_array[k]),op_depth_p[j][i],albedo))

    #Crea funcion fuente y temperatura en 3D
    funcion_S = RegularGridInterpolator((z_array, y_array, x_array), Sfunction,bounds_error=False,fill_value=None)
    funcion_T = RegularGridInterpolator((z_array, y_array, x_array), Temfunction,bounds_error=False,fill_value=None)
    return funcion_S, funcion_T

#--------------------------------------------------------------------------------
def FUN_f(t,tau,alb):
    eps = np.sqrt(1.-alb)
    fff = np.exp(-np.sqrt(3.)*eps*t) + np.exp(np.sqrt(3.)*eps*(t-tau))
    fff = fff/( np.exp(-np.sqrt(3.)*eps*tau)*(eps-1.) - (eps+1.) )
    return fff

#--------------------------------------------------------------------------------
#Lee las tablas de opacidad DSHARP
#Load opacities
with np.load('default_opacities_smooth.npz') as d:
    a_w     = d['a']
    gsca_w  = d['g']
    lam_w   = d['lam']
    k_abs_w = d['k_abs']
    k_sca_w = d['k_sca']

lam_avgs = wl
# We split the opacities within the range of frequency to make the calculations faster
k_abs_w = k_abs_w[(0.9*lam_avgs<lam_w) & (1.1*lam_avgs>lam_w),:]
k_sca_w = k_sca_w[(0.9*lam_avgs<lam_w) & (1.1*lam_avgs>lam_w),:]
k_sca_w = k_sca_w*(1. -  gsca_w[(0.9*lam_avgs<lam_w) & (1.1*lam_avgs>lam_w),:])
lam_w = lam_w[(0.9*lam_avgs<lam_w) & (1.1*lam_avgs>lam_w)]

opac_grid = opacity.size_average_opacity(lam_avgs, a_w, lam_w, k_abs_w.T, k_sca_w.T, q=3.5, plot=True)


function_ext = interpolate.interp1d(a_w, opac_grid['ka'][:]+opac_grid['ks'][:],kind='cubic')
function_alb = interpolate.interp1d(a_w, opac_grid['ks'][:]/(opac_grid['ka'][:]+opac_grid['ks'][:]),kind='cubic')
if not scattering:
    function_alb = interpolate.interp1d(a_w, np.zeros((np.shape(opac_grid['ks'][:]))),kind='cubic')
