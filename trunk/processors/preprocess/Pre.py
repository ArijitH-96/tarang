
######################################################################
######################################################################
####### @ Author : Shubhadeep Sadhukhan ##############################
########## Pre-Processing script TARANG ##############################
######################################################################
from pylab import*
import yaml
with open("para.yaml",'r') as stream1:
    data=(yaml.load(stream1))
kfactor1=data["field"]["kfactor"][0]
kfactor2=data["field"]["kfactor"][1]
kfactor3=data["field"]["kfactor"][2]

N1=data["field"]["N"][0]
N2=data["field"]["N"][1]
N3=data["field"]["N"][2]

L1=2*pi*data["field"]["kfactor"][0]
L2=2*pi*data["field"]["kfactor"][1]
L3=2*pi*data["field"]["kfactor"][2]

wavenumber_switch=data["field"]["waveno_switch"]
if N2>1:
    L=(L1*L2*L3)**(1/3.)
    N=(N1*N2*N3)**(1/3.)
else:
    L=(L1*L3)**(1/2.)
    N=(N1*N3)**(1/2.)
dissipation_coefficient=data["field"]["diss_coefficients"][0]
epsilon=1 # Anas: data["io"]["double_para"][0]
#N=int(N+0.01)

##########################################################################
##########################################################################
def Min_radius_outside():
    if wavenumber_switch:
        return int (ceil(sqrt( ((N1/2)*kfactor1)**2+ (N2/2*kfactor2)**2 + (N3/2*kfactor3)**2)));
    else:
        return int (ceil(sqrt((N1/2)**2 + (N2/2)**2 + (N3/2)**2 )));

def Number_of_shell():
    min_radius=Min_radius_outside()
    return min_radius+1
N=Number_of_shell()
print (N)
def print_as_blitz(ek,N,initial):
    N=N-1
    initial.write("(0,%d)\n"%N)
    for i in range(N+1):
        if i==0:
            initial.write("[ %e "%ek[i])
        elif i==N:
            initial.write("%e ]\n"%ek[i])
        else:
            initial.write("%e "%ek[i])
##########################################################################
def Pope_Model(dissipation_coefficient,epsilon):
    p0=2;
    cL=6.78;
    ceta=0.40;
    beta=5.2;
    Kolmogorov_constant=1.5;
    eta=0.42*(dissipation_coefficient**3/epsilon)**0.25;
    ek=zeros(N)
    for i in range(1,N-1):
        k=i
        fL_arg=k*L/sqrt((k*L/2)**2+cL**2)
        fL=fL_arg**(5./3+p0)
        feta_arg=((k*eta)**4+ceta**4)**0.25-ceta
        feta=exp(-beta*feta_arg)
        ek[i]=Kolmogorov_constant*epsilon**(2./3)*k**(-5./3)*fL*feta
    return ek



def Pao_Model(dissipation_coefficient,epsilon):
    Kolmogorov_constant=1.5;
    ek=zeros(N)
    eta=(dissipation_coefficient**3/epsilon)**0.25
    for i in range(1,N-1):
        k=i

        ek[i]=Kolmogorov_constant*epsilon**(2./3)*k**(-5/3.)*exp(-(k*eta)**(4/3.))
    return ek
##################################################
def energy_init():
    x=zeros(N)
    x=Pope_Model(dissipation_coefficient,epsilon)
    E=2
    x=E*x/sum(x)
    print ("Total energy",sum(x))
    return x
def helicity_init():
    x=zeros(N)
    Ek = energy_init()
    for i in range(N):
        x[i]= 0 #0.5*i*Ek[i]
    print ("Total helicity",sum(x))
    return x
def temperature_init():
    x=zeros(N)
    for i in range(N):
        x[i]=0.0
    print ("Total temperature energy",sum(x))
    return x
def magnetic_energy_init():
    x=zeros(N)  
    for i in range(N):
        if 1<=i<=3:
            x[i]=0.0001
    print ("Total magnetic energy",sum(x))
    return x
def magnetic_helicity_init():
    x=zeros(N)
    Hm = magnetic_energy_init()
    for i in range(N):
        x[i] = 0#0.5*i*Hm[i]
    print ("Total magnetic helicity",sum(x))
    return x
def magnetic_crosshelicity_init():
    x=zeros(N)
    for i in range(N):
        x[i]=0.0
    print ("Total cross helicity",sum(x))
    return x
####################################################

def energy_force():
    x=zeros(N)
    for i in range(N):
        if 3<=i<4:
            x[i]=0.0
    print ("Total energy supply",sum(x))
    return x
def helicity_force():
    x=zeros(N)
    khf = energy_force()
    for i in range(N):
        if 3<=i<=6:
            x[i]= 0 #0.5*i*khf[i]
    print ("Total helicity supply",sum(x))
    return x
def temperature_force():
    x=zeros(N)
    for i in range(N):
        x[i]=0.0
    print ("Total thermal supply",sum(x))
    return x
def magnetic_energy_force():
    x=zeros(N)
    for i in range(N):
        x[i]=0.0
    print ("Total magnetic energy supply",sum(x))
    return x
def magnetic_helicity_force():
    x=zeros(N)
    for i in range(N):
        x[i]=0.0
    print ("Total magnetic helicity supply",sum(x))
    return x
def magnetic_crosshelicity_force():
    x=zeros(N)
    for i in range(N):
        x[i]=0.0
    print ("Total cross helicity supply",sum(x))
    return x
#################################################
initial=open("initial.txt","w")
forcing=open("forcing.txt","w")
def ifluid():
    ek=energy_init()
    hk=helicity_init()
    eps=energy_force()
    epsh=helicity_force()
    '''
    loglog(ek,lw=3)
    xlabel("k",fontsize=15)
    ylabel("Ek",fontsize=15)
    title("Pao Spectrum",fontsize=15)
    savefig("Pao.jpg")

    show()
    '''
    print_as_blitz(ek,N,initial)
    print_as_blitz(hk,N,initial)
    print_as_blitz(eps,N,forcing)
    print_as_blitz(epsh,N,forcing)

def RBC():
    ek=energy_init()
    hk=helicity_init()
    T=temperature_init()
    eps=energy_force()
    epsh=helicity_force()
    epsT=temperature_force()
    print_as_blitz(ek,N,initial)
    print_as_blitz(hk,N,initial)
    print_as_blitz(T,N,initial)
    print_as_blitz(eps,N,forcing)
    print_as_blitz(epsh,N,forcing)
    print_as_blitz(epsT,N,forcing)
def imhd():
    ek=energy_init()
    hk=helicity_init()
    eb=magnetic_energy_init()
    hb=magnetic_helicity_init()
    hc=magnetic_crosshelicity_init()
    eps=energy_force()
    epsh=helicity_force()
    epsb=magnetic_energy_force()
    epshb=magnetic_helicity_force()
    epshc=magnetic_crosshelicity_force()
    print_as_blitz(ek,N,initial)
    print_as_blitz(hk,N,initial)
    print_as_blitz(eb,N,initial)
    print_as_blitz(hb,N,initial)
    print_as_blitz(hc,N,initial)
    print_as_blitz(eps,N,forcing)
    print_as_blitz(epsh,N,forcing)
    print_as_blitz(epsb,N,forcing)
    print_as_blitz(epshb,N,forcing)
    print_as_blitz(epshc,N,forcing)


program_kind=data["program"]["kind"]
if program_kind=="FLUID_INCOMPRESS":
    ifluid()
elif program_kind=="RBC":
    RBC()
elif program_kind=="MHD_INCOMPRESS":
    imhd()
else:
    print ("Not Implemented")
