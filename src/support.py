import numpy as np
import matplotlib.pylab as plt
# support routines for jupyter notebook

#astro distance
AU = 1.496*1e11 #m


#some constants
sigmaSB = 5.670373e-8 #W⋅m^-2⋅K^-4
h = 6.626e-34
c = 3.0e+8
k = 1.38e-23

Rsun = 6.957*1e8 # m
Tsun = 5778 #K
Lsun = 4. * np.pi * Rsun**2 *  sigmaSB * Tsun**4 # Watts

#distances of planets to sun

d_mercury = 0.39* AU #m
d_venus   = 0.723* AU #m
d_earth   = 1* AU #m
d_mars    = 1.524* AU #m
d_jupiter = 5.203* AU #m
d_saturn  = 9.539* AU #m
d_uranus  = 19.18* AU #m
d_neptune = 30.06* AU #m

pnames     = np.array(['Sun','Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune'])
pcolor     = np.array(['yellow','gray', 'lemonchiffon', 'blue', 'darkred', 'orange', 'palegoldenrod', 'royalblue', 'cornflowerblue'])
pdistances = np.array([0.,d_mercury, d_venus, d_earth, d_mars, d_jupiter, d_saturn, d_uranus, d_neptune])
pradii     = 0.5*np.array([Rsun*2./1000.,4879., 12104., 12756., 6792., 142984., 120536., 51118., 49528.])*1000.


def planck(wav, T):
    a = 2.0*h*c**2
    b = h*c/(wav*k*T)
    intensity = a/ ( (wav**5) * (np.exp(b) - 1.0) )
    return intensity

def plot_hz(lower_limit=1,upper_limit=2):
    if(upper_limit > 1e8):
        #assume results are given in meter
        upper_limit/=AU
        lower_limit/=AU
        
    x = pdistances/AU
    y = np.ones(len(x))

    fig, ax = plt.subplots(figsize=(14, 7))
    ax.set_title('The habitable zone in the solar system', fontsize=20)
    ax.set_ylim([0.5,1.5])
    ax.set_xlim([0.25,40])
    ax.set_xscale('log')
    ax.set_xlabel(r'$d$ [AU]', fontsize=20)
    ax.get_yaxis().set_visible(False)
    if(upper_limit>lower_limit):
        ax.axvspan(lower_limit, upper_limit, alpha=0.3, color='green')
    s0=400	
    for i, txt in enumerate(pnames):
        ax.scatter(x[i], y[i], color=pcolor[i], s=s0*pradii[i]/pradii[2], marker='o')    
        ax.annotate(txt, (x[i], y[i]+0.1), fontsize=10)

def plot_HR(L=10**np.random.uniform(-5,6,200)*Lsun,T=np.random.uniform(0.25,10,200)*Tsun,Tmin=1000,Tmax=40000,Lmin=1e-5,Lmax=1e6):
    y, x = np.mgrid[np.log10(Lmin):np.log10(Lmax):100j, np.log10(Tmin):np.log10(Tmax):100j]
    size = np.sqrt(10**(y)*(Tsun/(10**x))**4)
    temp = x
    fig, ax = plt.subplots(figsize=(10, 12))
    ax.set_ylim([y.min(),y.max()])
    ax.set_xlim([x.min(),x.max()])
    
    lgLsun=np.log10(1)
    lgTsun=np.log10(Tsun)
    
    lgL=np.log10(L/Lsun)
    lgT=np.log10(T)
    R=np.sqrt(10**(lgL)*(Tsun/(10**lgT))**4)
    
    
    levels = [0.001,0.01,0.1,1,10,100,1000]
    CS = ax.contour(x, y, size,colors='k',levels=levels)
    
    
    fmt = {}
    strs = [str(levels[0]) + r' $R_{\odot}$ ', str(levels[1]) + r' $R_{\odot}$ ', str(levels[2]) + r' $R_{\odot}$ ', str(levels[3]) + r' $R_{\odot}$ ', str(levels[4]) + r' $R_{\odot}$ ', str(levels[5]) + r' $R_{\odot}$ ', str(levels[6]) + r' $R_{\odot}$ ']
    for l, s in zip(CS.levels, strs):
        fmt[l] = s
    
    
    plt.clabel(CS, inline=1,fmt=fmt, fontsize=10)

    s0=50#size of sun
    
    ax.set_xlabel(r'surface temperature $T_{\star}$ [K]', fontsize=20)
    ax.set_ylabel(r'luminosity $L_{\star}$ [$L_{_\odot}$]', fontsize=20)
    ax.set_title('Hertzsprung–Russell diagram', fontsize=20)
    ax.scatter(lgT,lgL,s=R*s0,c=T,cmap='magma_r')
    #ax.scatter(lgTsun,lgLsun,c=T,s=s0,cmap='magma_r')
    ax.invert_xaxis()

    
def gen_spectra():
    
    Ts,Mag=np.loadtxt('./data/stars.txt', usecols=[0,1], unpack=True)
    Ls=Lsun*10**(0.4*(4.74-Mag))
    #find white dwarf
    foo =np.where((Ls<5e-2*Lsun)*(Ts>10**3.8))

    Twd=Ts[foo][0]
    Lwd=Ls[foo][0]
    
    #find red giant
    foo =np.where((Ls>1e2*Lsun)*(Ts<10**3.6))
    Trg=Ts[foo][0]
    Lrg=Ls[foo][0]
    
    # generate x-axis in increments from 1nm to 3 micrometer in 1 nm increments
    # starting at 1 nm to avoid wav = 0, which would result in division by zero.
    wavelengths = np.arange(1e-9, 2e-6, 1e-10)               
    
    
    sigma=0.5
                  
    #generate spectrum sun:
    spec_sun=planck(wavelengths, Tsun)
    
    mx=np.max(spec_sun)
    spec_sun+=mx*0.01*np.random.normal(0, 1, len(wavelengths))
    
    #generate spectrum white dwarf:
    spec_wd=planck(wavelengths, Twd)
    
    mx=np.max(spec_wd)
    spec_wd+=mx*0.01*np.random.normal(0, 1, len(wavelengths))          
    
    #generate spectrum red giant:
    spec_rg=planck(wavelengths, Trg)
    
    mx=np.max(spec_rg)
    spec_rg+=mx*0.01*np.random.normal(0, 1, len(wavelengths))          
    
    np.savez('./data/sun_data',wavelengths=wavelengths*1e9,spectrum=spec_sun,T=Tsun,L=Lsun)
    np.savez('./data/whitedwarf_data',wavelengths=wavelengths*1e9,spectrum=spec_wd,T=Twd,L=Lwd)
    np.savez('./data/redgiant_data',wavelengths=wavelengths*1e9,spectrum=spec_rg,T=Trg,L=Lrg)

    
def plot_data_spec(objectnr):
    
    if(objectnr==1):
        data=np.load('./data/sun_data.npz')
    if(objectnr==2):
        data=np.load('./data/whitedwarf_data.npz')
    if(objectnr==3):
        data=np.load('./data/redgiant_data.npz')

    x = data['wavelengths']
    y = data['spectrum']

    plt.figure(figsize=(14, 7))
    plt.plot(x,y,lw=0.5)
    plt.xlabel(r'wavelength $\lambda$ [nm]')
    plt.ylabel(r'SSI [W $m^{-2}\, nm^{-1}$]')
    plt.grid()
    plt.show()


def get_spectra():

    wavelengths = np.arange(1e-9, 2e-6, 1e-10) 

    Temps=np.array([7000,6500,6000,5500,5000,4500,4000])
    Temp=[]
    intensity=[]
    wavelength=[]
    plt.figure(figsize=(14, 7))
    for T in Temps:
        B = planck(wavelengths, T)
        intensity.append(B)
        Temp.append(T)
        wavelength.append(wavelengths*1e9)
        plt.plot(wavelengths*1e9, B,label='T='+'{:06.2f}'.format(T)+' [K]') 
        plt.ylabel(r'$B(\lambda)$')
        plt.xlabel(r'$\lambda$')

    # show the plot
    plt.grid()
    plt.legend()
    plt.show()
    
    return np.array(wavelength),np.array(intensity),np.array(Temp)
    
    
    
    
    




