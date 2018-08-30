import numpy as np
import matplotlib.pylab as plt
# support routines for jupyter notebook

#astro distance
AU = 1.496*1e11 #m


#some constants
sigmaSB = 5.670373e-8 #W⋅m^-2⋅K^-4

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





