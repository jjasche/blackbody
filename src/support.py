import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
import matplotlib.ticker as ticker
# support routines for jupyter notebook

#astro distance
AU = 1.496*1e11 #m


#some constants
sigmaSB = 5.670373e-8 #W m^-2 K^-4
h = 6.626e-34
c = 3.0e+8
k = 1.38e-23
b = 2.8977729*1e-3# m⋅K


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

pnames_random =  np.array(['GC-818','Proxima Norma', 'Aporia', 'Nyota', 'Cantor', 'Algirae', 'Tarvos', 'Yuma', 'Gratian'])   
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
    ax.set_title('The habitable zone around a star', fontsize=20)
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

def plot_HR(Ls=10**np.random.uniform(-5,6,200)*Lsun,Ts=np.random.uniform(0.25,10,200)*Tsun,Tmin=1000,Tmax=40000,Lmin=1e-5,Lmax=1e6):
    
    T,Mag=np.loadtxt('./data/stars.txt', usecols=[0,1], unpack=True)
    L=Lsun*10**(0.4*(4.74-Mag))
    
    
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
    
    
    lgLs=np.log10(Ls/Lsun)
    lgTs=np.log10(Ts)
    Rs=np.sqrt(10**(lgLs)*(Tsun/(10**lgTs))**4)
    
    
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
    ax.scatter(lgTs[0],lgLs[0],s=s0*4,color='green')
    ax.scatter(lgTs[0],lgLs[0],s=s0*4,color='goldenrod',marker=(5, 2))
    ax.annotate("Object Nr. 1", (lgTs[0]+.1, lgLs[0]+.3), fontsize=18,color='green')
    ax.scatter(lgTs[1],lgLs[1],s=s0*4,color='red')
    ax.scatter(lgTs[1],lgLs[1],s=s0*4,color='goldenrod',marker=(5, 2))
    ax.annotate("Object Nr. 2", (lgTs[1]+.1, lgLs[1]+.3), fontsize=18,color='red')
    ax.scatter(lgTs[2],lgLs[2],s=s0*4,color='blue')
    ax.scatter(lgTs[2],lgLs[2],s=s0*4,color='goldenrod',marker=(5, 2))
    ax.annotate("Object Nr. 3", (lgTs[2]+.1, lgLs[2]+.3), fontsize=18,color='blue')
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

    
    fig =plt.figure(figsize=(14, 7)) 
    ax = fig.add_subplot(111)
  
    ax.plot(x,y,lw=0.5)
    ax.set_xlim([0.,1500])
    ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(25))
    ax.set_xlabel(r'wavelength $\lambda$ [nm]')
    ax.set_ylabel(r'SSI [W $m^{-2}\, nm^{-1}$]')
    ax.grid(b=True, which='major', linestyle='-')
    ax.grid(b=True, which='minor', linestyle=':')
    plt.show()
    
    return x, y

def get_object_luminosities():
    
    data=np.load('./data/sun_data.npz')
    L1=data['L']
    data=np.load('./data/whitedwarf_data.npz')
    L2=data['L']
    data=np.load('./data/redgiant_data.npz')
    L3=data['L']
    
    return L1, L2, L3

def get_object_luminosities_and_temperatures():
    
    data=np.load('./data/sun_data.npz')
    L1=data['L']
    T1=data['T']
    data=np.load('./data/whitedwarf_data.npz')
    L2=data['L']
    T2=data['T']
    data=np.load('./data/redgiant_data.npz')
    L3=data['L']
    T3=data['T']
    
    return T1,T2,T3,L1, L2, L3






def get_spectra():

    wavelengths = np.arange(1e-10, 2e-6, 0.1e-10) 

    Temps=np.array([8000,7500,7000,6500,6000,5500,5000,4500,4000,3000])
    Temp=[]
    intensity=[]
    wavelength=[]
    
    fig =plt.figure(figsize=(14, 7)) 
    ax = fig.add_subplot(111)
  
    for T in Temps:
        B = planck(wavelengths, T)
        intensity.append(B)
        Temp.append(T)
        wavelength.append(wavelengths*1e9)
        ax.plot(wavelengths*1e9, B,label='T='+'{:06.2f}'.format(T)+' [K]') 
        #ax.xticks(np.arange(np.min(wavelengths*1e9), np.max(wavelengths*1e9), step=100))
        ax.set_title('Black Body spectra')
        ax.set_ylabel(r'SSI [W $m^{-2}\, nm^{-1}$]')
        ax.set_xlabel(r'wavelength $\lambda [\mathrm{nm}]$')

    # show the plot
    ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(25))
    ax.set_xlim(0,1500)
    ax.grid(b=True, which='major', linestyle='-')
    ax.grid(b=True, which='minor', linestyle=':')
    ax.legend()
    plt.show()
    
    return np.array(wavelength),np.array(intensity),np.array(Temp)



# calculate temperature of planet
def planet_temperature(Lstar=Lsun):
    
    Rp = np.linspace(0,6*AU, 200)
    Tp = (0.25*Lstar/(4*np.pi*sigmaSB*Rp**2))**(0.25)
    
   
    
    fig =plt.figure(figsize=(14, 7)) 
    ax = fig.add_subplot(111)
    ax.set_xlim([0,6.])
    ax.set_ylim([0,800.])
    ax.plot(Rp/AU,Tp,color='blue')
    ax.set_ylabel(r'$T \, [\mathrm{K}]$')
    ax.set_xlabel(r'$R_p \, [\mathrm{AU}]$')
    ax.grid(b=True, which='major', linestyle='-')
    ax.grid(b=True, which='minor', linestyle=':')
    ax.set_title(r'Temperature of planets at distance  $R_p$')
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.25))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.125))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(50))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(25))
    
    plt.show()
    

    
#SOLUTIONS

#1) Solution Wien's displacement
def solution_wien_displacement(lambdas=None):
    Tmax=np.array([8000,7500,7000,6500,6000,5500,5000,4500,4000,3000])
    lmax = b/Tmax*1e9 #nm
    
    plt.figure(figsize=(14, 7))    
    plt.ylim([0,10000])
    plt.xlim([350,800])
    plt.scatter(np.array(lmax),np.array(Tmax))
    plt.plot(np.array(lmax),np.array(Tmax))
    plt.ylabel(r'$T \, [\mathrm{K}]$')
    plt.xlabel(r'$\lambda \, [\mathrm{nm}]$')
    
    plt.title('Wien\'s displacement law')
    Temps=None
    Lums=None
    if lambdas is not None:
        Temps = b*1e9/lambdas #K
        colors = cm.jet(np.linspace(0, 1, len(lambdas)))
        for l in np.arange(len(lambdas)):
            plt.scatter(lambdas[l],Temps[l],color=colors[l])
            plt.plot([lambdas[l],lambdas[l]],[0,Temps[l]],'--',color=colors[l])
            plt.plot([0,lambdas[l]],[Temps[l],Temps[l]],'--',color=colors[l])
            plt.annotate("T = {0:.2f} [K]".format(Temps[l]), (lambdas[l]+5., Temps[l]+10.), fontsize=18,color=colors[l])
    
        #check Solution
        T1, T2, T3, L1, L2, L3 = get_object_luminosities_and_temperatures()
        
        print (T2,Temps[1])
        
        dTemp=200 #K  read-of error

        if np.fabs(Temps[0]-T1)<dTemp:
            print ('Temperature of Stellar object Nr. 1 identified correctly')
            print (r'Total Luminosity of Stellar object Nr. 1 $L_{\star}$=',L1/Lsun)
        else:
            print ('Temperature of Stellar object Nr. 1 does not match observations')
            print ('Try again')
            
        if np.fabs(Temps[1]-T2)<dTemp:
            print ('Temperature of Stellar object Nr. 2 identified correctly')
            print (r'Total Luminosity of Stellar object Nr. 2 $L_{\star}$=',L2/Lsun)
        else:
            print ('Temperature of Stellar object Nr. 2 does not match observations')
            print ('Try again')
            
        if np.fabs(Temps[2]-T3)<dTemp:
            print ('Temperature of Stellar object Nr. 3 identified correctly')
            print (r'Total Luminosity of Stellar object Nr. 3 $L_{\star}$=',L3/Lsun)
        else:
            print ('Temperature of Stellar object Nr. 3 does not match observations')
            print ('Try again')
            
        Lums=np.array([L1,L2,L3])    
    
    plt.show()
    
    
    
    return Temps, Lums
    

#1) Solution Wien's displacement
def plot_solution_wien_displacement():    
    
    data=np.load('input_data_wien_displacement.npz')
    T_list=data['T_list']
    lambda_list=data['lambda_list']
    
    fig =plt.figure(figsize=(14, 7)) 
    ax = fig.add_subplot(111)
    
    ax.set_ylabel('Wavelength (nm)')
    ax.set_xlabel('Temperature (K)')
    ax.grid(b=True, which='major', linestyle='-')
    ax.grid(b=True, which='minor', linestyle=':')
    ax.set_title(r"Wien's displacement law ")
    
    Tmax=np.linspace(3000,8000,100)
    lmax = b/Tmax*1e9 #nm
    ax.set_xlim(3000,8000)
    ax.set_ylim(0,1000)
    ax.plot(T_list,lambda_list, linestyle=' ', marker='o', color='r', label='Square')
    ax.plot(Tmax,lmax, color='blue')
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(200))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(100))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(50))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(25))
    
    plt.show






