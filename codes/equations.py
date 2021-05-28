#!/usr/bin/env python

def EQ_density(x_d,y_d,z_d): #Dust density
    r = np.sqrt(x_d**2+y_d**2)
    H = EQ_Height(r)    
    if (r < Rin or r> Rout or abs(z_d)>2.5*H):                    #at z=2H the disk mass is 95.45%
        return 0.
    else:
        Sigma = 1.*r**(-1.)*np.exp(-r/Rc)
        return Sigma*np.exp(-0.5*(z_d/H)**2)/np.sqrt(2.*np.pi)/H #AU is not included due to integration is in terms of H


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


def EQ_amax(x,y,z):
    return 0.1 #In cm units

