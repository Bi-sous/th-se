import matplotlib.pyplot as plt
import numpy as np

lamb = 30

def bohren_battan_formula(n_i, n_m, F):
  e_m = n_m**2
  e_i = n_i**2
  b = 2*e_m/(e_i - e_m)*(e_i/(e_i-e_m)*np.log(e_i/e_m) - 1)
  return np.sqrt(((1-F)*e_m + F*b*e_i)/(1 - F + b*F))

def permice():
    t = 273-273. # celcius
    l = lamb/10. # cm
    sig0 = 18.8496e10
    einf = 3.168
    a = 0.288 + 0.0052*t + 0.00023*t**2
    sig = 1.26*np.exp(-12500.0/(1.9869*(t+273.)))
    ls = 9.990288e-5*np.exp(13200.0/(1.9869*(t+273)))
    #ls = 9.990288e-4*np.exp(13200.0/(1.9869*(t+273))) --> e-4 --> e-5 error pointed by ogushi
    rl = ls/l
    es = 203.168 + 2.5*t + 0.15*t**2
    denom = 1. + 2*rl**(1-a)*np.sin(a*np.pi/2.) + rl**(2*(1-a))
    ep = einf + (es - einf)*(1 + rl**(1-a)*np.sin(a*np.pi/2))/denom
    epp = (es - einf)*rl**(1-a)*np.cos(a*np.pi/2)/denom + sig*l/sig0
    return np.sqrt(ep + 1.0j*epp)

def permwat():
    t = 273-273. # celcius
    l = np.atleast_1d(lamb/10.) # cm
    sig = 12.5664e8
    sig0 = 18.8496e10
    einf = 5.27137 + 0.0216474*t - 0.00131198*t**2
    a = -16.8129/(t+273) + 0.0609265
    ls = 0.00033836*np.exp(2513.98/(t+273.))
    rl = ls/l
    es = 78.54*(1. - 4.579e-3*(t-25) + 1.19e-5*(t-25)**2 - 2.8e-8*(t-25)**3)
    denom = 1. + 2*rl**(1-a)*np.sin(a*np.pi/2.) + rl**(2*(1-a))
    ep = einf + (es - einf)*(1 + rl**(1-a)*np.sin(a*np.pi/2))/denom
    epp = (es - einf)*rl**(1-a)*np.cos(a*np.pi/2)/denom + sig*l/sig0
    plus = np.zeros_like(15)    
    epp += plus
    return  np.sqrt(ep + 1.0j*epp)


n_ice = permice()
n_wat = permwat()
n_air = np.sqrt(1.0005364)

dens_w = 1000
dens_i = 931
dens_a = 1.292

density = 100
vol = 4/3 * np.pi * (5e-3)**3
mass = density * vol

F_i_s = (density-dens_a)/(dens_i-dens_a)*100
F_a_s = 100-F_i_s

F_w = np.arange(0,100.1,0.1)/100
dens_f = np.linspace(100,1000, len(F_w))
F_i = (dens_f-dens_a+F_w*dens_a-F_w*dens_w)/(dens_i-dens_a)
F_a = 1-F_i-F_w


# type (a)
wat_in_ice = bohren_battan_formula(n_wat, n_ice, F_w/(F_w+F_i))
wat_in_ice_in_air = bohren_battan_formula(wat_in_ice, n_air, (F_w+F_i)/(F_w+F_i+F_a))
wat_in_ice_in_air_p = bohren_battan_formula(n_air, wat_in_ice, (F_a)/(F_w+F_i+F_a))

# type (b)
ice_in_wat = bohren_battan_formula(n_ice, n_wat, F_i/(F_w+F_i))
ice_in_wat_in_air = bohren_battan_formula(ice_in_wat, n_air, (F_w+F_i)/(F_w+F_i+F_a))
ice_in_wat_in_air_p = bohren_battan_formula(n_air, ice_in_wat, (F_a)/(F_w+F_i+F_a))

# type (c)
air_in_ice = bohren_battan_formula(n_air, n_ice, F_a/(F_a+F_i))
air_in_ice_in_wat = bohren_battan_formula(air_in_ice, n_wat, (F_a+F_i)/(F_w+F_i+F_a))
air_in_ice_in_wat_p = bohren_battan_formula(n_wat, air_in_ice, (F_w)/(F_w+F_i+F_a))

# type (d)
ice_in_air = bohren_battan_formula(n_ice, n_air, F_i/(F_i+F_a))
ice_in_air_in_wat = bohren_battan_formula(ice_in_air, n_wat, (F_i+F_a)/(F_i+F_a+F_w))
ice_in_air_in_wat_p = bohren_battan_formula(n_wat, ice_in_air, (F_w)/(F_i+F_a+F_w))

# type (e)
air_in_wat = bohren_battan_formula(n_air, n_wat, F_a/(F_a+F_w))
air_in_wat_in_ice = bohren_battan_formula(air_in_wat, n_ice, (F_a+F_w)/(F_a+F_w+F_i))
air_in_wat_in_ice_p = bohren_battan_formula(n_ice, air_in_wat, (F_i)/(F_a+F_w+F_i))

# type (f)
wat_in_air = bohren_battan_formula(n_wat, n_air, F_w/(F_w+F_a))
wat_in_air_in_ice = bohren_battan_formula(wat_in_air, n_ice, (F_w+F_a)/(F_a+F_w+F_i))
wat_in_air_in_ice_p = bohren_battan_formula(n_ice, wat_in_air, (F_i)/(F_a+F_w+F_i))

perm_drysnow = bohren_battan_formula(n_ice, n_air, F_i_s/100)
perm_drysnow = np.ones_like(dens_f)*perm_drysnow

perm_wat = np.ones_like(dens_f)*n_wat


#graph
plt.figure(dpi=200)

plt.subplot(221)
plt.plot(dens_f, np.real(air_in_wat_in_ice**2), 'b', label='(e)')
plt.plot(dens_f, np.real(air_in_wat_in_ice_p**2), 'r', label='(e )')
plt.grid()
plt.legend()
plt.xlim(100,1000)
plt.ylim(0,45)
plt.plot(dens_f, np.real(perm_drysnow**2), 'k--')
plt.plot(dens_f, np.real(perm_wat**2), 'k-.')
plt.subplot(222)
plt.plot(dens_f, np.imag(air_in_wat_in_ice**2), 'b', label='(e)')
plt.plot(dens_f, np.imag(air_in_wat_in_ice_p**2), 'r', label='(e )')
plt.legend()
plt.xlim(100,1000)
plt.ylim(0,45)
plt.plot(dens_f, np.imag(perm_drysnow**2), 'k--')
plt.plot(dens_f, np.imag(perm_wat**2), 'k-.')
plt.grid()

plt.subplot(223)
plt.plot(dens_f, np.real(wat_in_air_in_ice**2), 'b', label='(f)')
plt.plot(dens_f, np.real(wat_in_air_in_ice_p**2), 'r', label='(f )')
plt.grid()
plt.legend()
plt.xlim(100,1000)
plt.ylim(0,45)
plt.plot(dens_f, np.real(perm_drysnow**2), 'k--')
plt.plot(dens_f, np.real(perm_wat**2), 'k-.')
plt.subplot(224)
plt.plot(dens_f, np.imag(wat_in_air_in_ice**2), 'b', label='(f)')
plt.plot(dens_f, np.imag(wat_in_air_in_ice_p**2), 'r', label='(f )')
plt.legend()
plt.xlim(100,1000)
plt.ylim(0,45)
plt.plot(dens_f, np.imag(perm_drysnow**2), 'k--')
plt.plot(dens_f, np.imag(perm_wat**2), 'k-.')
plt.grid()

#
#plt.plot(dens_f, np.real(ice_in_wat_in_air**2), 'b--', label='(b)')
#plt.plot(dens_f, np.real(air_in_ice_in_wat**2), 'r', label='(c)')
#plt.plot(dens_f, np.real(ice_in_air_in_wat**2), 'r--', label='(d)')
#plt.plot(dens_f, np.real(air_in_wat_in_ice**2), 'r-.', label='(e)')
#plt.plot(dens_f, np.real(wat_in_air_in_ice**2), 'b-.',label='(f)')
#
#plt.plot(dens_f, np.real(wat_in_ice_in_air_p**2),'r', alpha=0.5, label='(a2)')
#plt.plot(dens_f, np.real(ice_in_wat_in_air_p**2),'r--', alpha=0.5, label='(b2)')
#plt.plot(dens_f, np.real(air_in_ice_in_wat_p**2), label='(c2)')
#plt.plot(dens_f, np.real(ice_in_air_in_wat_p**2), label='(d2)')
#plt.plot(dens_f, np.real(air_in_wat_in_ice_p**2), label='(e2)')
#plt.plot(dens_f, np.real(wat_in_air_in_ice_p**2), label='(f2)')


#plt.title('Partie rééle de la permittivité')
#
#plt.subplot(212)
#plt.plot(dens_f, np.imag(wat_in_ice_in_air**2), label='(a)')
#plt.plot(dens_f, np.imag(ice_in_wat_in_air**2), label='(b)')
#plt.plot(dens_f, np.imag(air_in_ice_in_wat**2), label='(c)')
#plt.plot(dens_f, np.imag(ice_in_air_in_wat**2), label='(d)')
#plt.plot(dens_f, np.imag(air_in_wat_in_ice**2), label='(e)')
#plt.plot(dens_f, np.imag(wat_in_air_in_ice**2), label='(f)')
#plt.xlim(100,1000)
#plt.ylim(0,45)
#plt.plot(dens_f, np.imag(perm_drysnow**2), 'k--')
#plt.plot(dens_f, np.imag(perm_wat**2), 'k--')
#plt.grid()
##plt.title('Partie imaginaire de la permittivité')
