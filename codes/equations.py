#!/usr/bin/env python

#def EQ_density(x_d,y_d,z_d): #Dust density
#    r = np.sqrt(x_d**2+y_d**2)
#    H = EQ_Height(r)    
#    if (r < Rin or r> Rout or abs(z_d)>2.5*H):                    #at z=2H the disk mass is 95.45%
#        return 0.
#    else:
#        Sigma = 1.*r**(-1.)*np.exp(-r/Rc)  # change 1
#        return Sigma*np.exp(-0.5*(z_d/H)**2)/np.sqrt(2.*np.pi)/H #AU is not included due to integration is in terms of H


rho_p = 1.25*(u.g*(u.cm)**(-3)) 
#s_1mm = 10**(-3)*u.m
s_10 = 10*10**(-6)*u.m 
#s_100 = 100*10**(-6)*u.m
M = const.M_sun
alpha = 10**(-3)*u.dimensionless_unscaled 
rho_0 = 2.83*10**(-10)*(u.g*(u.cm)**(-3))
h_0 = 3.33*10**(-2)*u.AU 
p = -2.25*u.dimensionless_unscaled 
q = -0.5*u.dimensionless_unscaled 

"""
Function for dust density equation (31) Takeuchi_2002
"""
def h_g(r):
    """
    Height
    """
    return h_0*r**((q+3)/2)

def omega_KMid(r):
    """
    Keplerian angular velocity
    """
    return (const.G*M/((r*u.AU)**3))**0.5


def v_KMid(r):
    """
    Keplerian tangential velocity
    """
    return (r*u.AU)*omega_KMid(r)


def rho_g(r, z):
    """
    Gas density
    """
    return rho_0*(r**p)*np.exp(-z**2/(2*h_g(r)**2))


def c(r):
    """
    Sound velocity
    """
    return h_g(r)*omega_KMid(r)


def v_t(r):
    """
    Mean thermal velocity
    """
    return (8/np.pi)**0.5*c(r)


def T_s(r, z, s=s_10):
    """
    Stopping time for particle of size 10 µm
    """
    return (rho_p*s*v_KMid(r))/(rho_g(r, z)*(r*u.AU)*v_t(r))


def rho_d0(r):
    """
    Dust density in the mid plane in g*cm-2
    """
    return 0.01  # constant

"""
New density equation
"""
def EQ_density(x_d, y_d, z_d):
    """
    Dust density for all the disk. Equation (31) Takeuchi_2002
    """
    r = np.sqrt(x_d**2 + y_d**2)
    z_d = z_d*u.AU
    H = h_g(r)    
    if (r < Rin) or (r > Rout) or (abs(z_d) > 2.5*H):
        return 0.
    else:
        Sc = 1*u.dimensionless_unscaled  # particle well coupled to the gas
        fraction = 0.5*(z_d/H)**2
        a = -fraction  # first term inside exp
        # second term inside exp
        b = -Sc*T_s(r, 0*u.AU)*(np.exp(fraction)-1*u.dimensionless_unscaled)/alpha
        rho_d = rho_d0(r)*np.exp(a + b)
        return rho_d.value


def EQ_temperature(x_d,y_d,z_d):
    return funcion_T([z_d,y_d,x_d])
    #r = np.sqrt(x_d**2+y_d**2)
    #H = EQ_Height(r)
    #Temp = 10. (r/1.)**(-0.5)


def EQ_Height(r):
    if (r > Rin):
        return 0.1*r**(1.0) #In AU units
    else:
        return 0.1*(Rin)**(1.0) #In AU units


def f_modelo_gaussianas(radio, A, sigma):
    """
    Modelo de las emisiones que consiste en dos gaussianas de mismo ancho
    (mismo sigma)

    Input:
    ============
    params : [tuple] tupla con los parámetros de la amplitud y centro de las
    dos gaussianas, y el sigma que es el mismo para ambas
    x : [float] valor a evaluar
    ============

    Output:
    ============
    [float] valor de la curva en el punto x
    ============
    """
    #A1, mu1, sigma = 0.7, 50, 2
    mu1 = 50
    evaluacion = A * np.exp(-(radio-mu1)**2/sigma**2) + 0.1
    return evaluacion


def EQ_amax(x, y, z):
    radius = np.sqrt(x**2 + y**2)
    grain_size = f_modelo_gaussianas(radius, 0.7, 2)-\
                 (f_modelo_gaussianas(radius, 0.3, 5)-0.2)
    return 1  # In cm units
