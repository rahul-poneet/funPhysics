import matplotlib.pyplot as plt
import numpy as np
# ====================================================================================================
#
#
#                .'|''-.
#               |  |    '.
#              |   |      '.
#             |    |        '.        The tail looks wrong, but I can't figure out how to lift it
#            '     |          '.     /
#          .'      |            '.  /
#      __.'________|______________'-.___________
#
def maxwellBoltzmannWithVelocitySquared(fig,ax, v, molarMass,temperature,crayon='None'):
    # Phys3070 - Atmospheric escape Maxwell Boltzmann Distribution with v^2.
    # Parameter: v                         - Domain.
    # Parameter: molarMass                 - The molar mass of the kinetic molecule.
    # Parameter: temperature               - Temperature of the distribution.
    # Parameter: crayon                    - Titan eat crayon?
    # Returns:   dontDieBoltzmann          - A nparray of mapped y-axis heights for domain v.
    # Local variable: kB                   - Boltzmann's constant.
    # Local variable: T                    - Renamed parameter temperature.
    # Local variable: N_AVOGADRO           - Avagadro's Number.
    # Local variable: molarMassInKilograms - Kilogram molar mass of parameter molarMass.
    # Local variable: m                    - Kilogram mass of each molecule.
    # Local variable: A                    - Sigma component of amplitude.
    # Local variable: B                    - Pi component of amplitude.
    # Local variable: C                    - Velocity squared component of amplitude.
    # Local variable: numerator            - The numerator of the exponent from definition.
    # Local variable: denominator          - The denominator of the expoent from definintion.
    # Local variable: dontDieBoltzmann     - A nparray of mapped y-axis heights for domain v.
    # Local variable: molecule             - 
    # Local variable: xlabel               - 
    # Local variable: ylabel               - 
    # Local variable: title                - 
    kB = 1.38e-23# J/K
    T = temperature
    # mass
    N_AVOGADRO = 6.02214076e23# molecules/mole
    molarMassInKilograms = molarMass *1e-3# kg/mole
    m = molarMassInKilograms / N_AVOGADRO# kg/molecule
    # amplitude components
    A = ( m / (2*kB*T) )**(3/2)
    B = (4) / (np.pi)**(1/2)
    C = (v)**(2)
    # exp    
    numerator = m*( (v)**(2) )
    denominator = 2*(kB*T)
    dontDieBoltzmann= A*B*C* np.exp((-1)* numerator/denominator)

    molecule = ""
    if (molarMass == 2*15.9994):
        molecule = "Oxygen Gas (O2)"
    if (molarMass == 15.9994):
        molecule = "Oxygen Atom (O)"
    if (molarMass == 2*1.00794):
        molecule = "Hydrogen Gas (H2)"
    if (molarMass == 1.00794):
        molecule = "Hydrogen Atom (H)"
    label = molecule + " at " + str(T) + "K"
    #plots
    if (crayon == 'None'):
        ax.plot(v,dontDieBoltzmann,label=label)
    else:
        ax.plot(v,dontDieBoltzmann,c=crayon)
        ax.fill(v,dontDieBoltzmann,alpha=0.1,c=crayon,label=label)
    
    xlabel = "Velocity [m/s]"
    ylabel = "Probability"
    title = "Maxwell-Boltzmann Distribution with Squared Velocity\n\n F(v) = A*B*C* exp[ -mv^2 / 2kT ]\n\n A=(m/2kT)^(3/2);     B=4*(π)^-1/2;     C= v^2\n"
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_yticks([])
    ax.set_xlim([0,max(v)])
    ax.set_ylim([0,1.1*max(dontDieBoltzmann)])
    fig.set_figwidth(20)
    fig.set_figheight(7)
    fig.legend()

    return(dontDieBoltzmann)
# ====================================================================================================
#
#           . - ' ' - .   
#         ' _.-'''''-._ ' 
#       ' .'           '. '
#      ' /               \ '
#     ' |        _________| '      ----------- o         E S C A P E
#     ' |        R_EARTH  | '   ----- o                V E L O C I T Y
#     '  \               /  '       --------------- o
#      '  '.           .'  '
#        '  '-._____.-'  '
#           ' .  _  . '
#
def escapeVelocity():
    # Calculates the escape velocity of Earth at the surface (at R_EARTH) and returns it.
    # Parameters:     None
    # Returns:        v_escape    - The calculated escape velocity at R_EARTH.
    # Local Variable: R_EARTH     - The radius of the Earth.
    # Local Variable: M_EARTH     - The mass of the Earth.
    # Local Variable: G           - Gravitational constant.
    # Local Variable: numerator   - Numerator of the escape velocity.
    # Local Variable: denominator - Denominator of the escape velocity.
    # Local Variable: v_escape    - The calculated escape velocity at R_EARTH.
    R_EARTH = 0.64e7
    M_EARTH = 5.97e24# kg
    G = 6.67e-11  
    numerator = 2* G* M_EARTH
    denominator = R_EARTH
    v_escape = (numerator/denominator)**(1/2) # sqrt( 2GM/R )
    print("----- DEBUG OUTPUT from function escapeVelocity() ----- returning ",v_escape,".")
    return(v_escape)
# ====================================================================================================
#
#         .''-._   _.-''.
#        /  _--_\ / _--_ \
#        | |    \|/    | | 
#         \ \   |||   / / 
#          \_.-'''''-._/  Crustal Magnetosphere (check scale)
#         .' ''-..--'' '.  
#        /               \  
#       |        _________|        ----------- o         E S C A P E
#       |        R_MARS   |     ----- o                V E L O C I T Y
#        \               /          --------------- o
#         '. __......__.'   
#          /'-._____.-'\
#         / /   |||   \ \
#        | -_  _/|\_  _- |
#        \   '' / \ ''   / 
#         -__.-'   '-.__-
#
def escapeVelocityOfMars():
    # Calculates the escape velocity of Earth at the surface (at R_MARS) and returns it.
    # Parameters:     None
    # Returns:        v_escape    - The calculated escape velocity at R_MARS.
    # Local Variable: R_EARTH     - The radius of the Earth.
    # Local Variable: M_EARTH     - The mass of the Earth.
    # Local Variable: R_MARS      - The radius of Mars.
    # Local Variable: M_MARS      - The mass of the Mars.
    # Local Variable: G           - Gravitational constant.
    # Local Variable: numerator   - Numerator of the escape velocity.
    # Local Variable: denominator - Denominator of the escape velocity.
    # Local Variable: v_escape    - The calculated escape velocity at R_MARS.
    R_EARTH = 0.64e7
    M_EARTH = 5.97e24# kg
    R_MARS = 0.533 * R_EARTH
    M_MARS = 0.107 * M_EARTH# kg
    G = 6.67e-11
    numerator = 2* G* M_MARS
    denominator = R_MARS
    v_escape = (numerator/denominator)**(1/2) # sqrt( 2GM/R )
    print("----- DEBUG OUTPUT from function escapeVelocity() ----- returning ",v_escape,".")
    return(v_escape)
# ====================================================================================================
#
#
#
#
#
#
#
#
#
def plotEscapeVelocity(fig,ax,v_escape,yMax,crayon,label):
    ax.plot([v_escape,v_escape],[0,yMax],crayon,label=label)
    fig.legend()
    return()
def adjustPlotHeight(ax,yMax):
    ax.set_ylim([0,1.1*yMax])
    return()

# ====================================================================================================
# Main
# ====================================================================================================
# Crayons
crayon = ['#ff0000','#ff0055','#5500ff','#0000ff']# Titan eat crayon?
crayonForMolecularO2 = '#0000ff'#   blue   = O2
crayonForAtomicO = '#00bbbb'#       cyan   = O
crayonForMolecularH2 = '#770077'#   purple = H2
crayonForAtomicH = '#ff0000'#       red    = H


# ====================================================================================================
# O2 at 300K
# ====================================================================================================
T = 300# K
molarMass = 2* 15.9994# O2 [g/mole]
v_max = 1500# K
v = np.arange(0,v_max,1)
fig1,ax1=plt.subplots()
y1 = maxwellBoltzmannWithVelocitySquared(fig1,ax1,v,molarMass,T,crayonForMolecularO2)
# ====================================================================================================
# Include Atomic O
molarMass = 15.9994# O [g/mole]
y1 = np.append(y1, maxwellBoltzmannWithVelocitySquared(fig1,ax1,v,molarMass,T,crayonForAtomicO))
adjustPlotHeight(ax1,max(y1))
# ====================================================================================================


# ====================================================================================================
# H2 at 100 and 300 K
# ====================================================================================================
T = [100,300]
molarMass = 2* 1.00794# H2 [g/mole]
v_max = 8000# K
v = np.arange(0,v_max,1)
y2 = np.array(0)
i = 0
fig2,ax2=plt.subplots()
while(i<len(T)):
    y2 = np.append(y2, maxwellBoltzmannWithVelocitySquared(fig2,ax2,v, molarMass,T[i],crayonForMolecularH2))
    i+=1
# ====================================================================================================
# Include Atomic H
molarMass = 1.00794# H [g/mole]
i=0
while(i<len(T)):
    y2 = np.append(y2, maxwellBoltzmannWithVelocitySquared(fig2,ax2,v, molarMass,T[i],crayonForAtomicH))
    i+=1
adjustPlotHeight(ax2,max(y2))
# ====================================================================================================


# ====================================================================================================
# O2 and H2 at 300K
# ====================================================================================================
T = 300# K
molarMassO2 = 2* 15.9994# O2 [g/mole]
molarMassH2 = 2* 1.00794# H2 [g/mole]
v_max = 12000# K
v = np.arange(0,v_max,1)
fig3,ax3=plt.subplots()
y3 = maxwellBoltzmannWithVelocitySquared(fig3,ax3,v,molarMassO2,T,crayonForMolecularO2)
y3 = np.append(y3, maxwellBoltzmannWithVelocitySquared(fig3,ax3,v,molarMassH2,T,crayonForMolecularH2))
# ====================================================================================================
# Include Atoms
molarMassAtomicO = 15.9994# O [g/mole]
molarMassAtomicH = 1.00794# H [g/mole]
y3 = np.append(y3, maxwellBoltzmannWithVelocitySquared(fig3,ax3,v,molarMassAtomicH,T,crayonForAtomicH))
y3 = np.append(y3, maxwellBoltzmannWithVelocitySquared(fig3,ax3,v,molarMassAtomicO,T,crayonForAtomicO))
adjustPlotHeight(ax3,max(y3))
# ====================================================================================================


# ====================================================================================================
# H2 at 25, 100, 300, 1200 K
# ====================================================================================================
T = [25,100,300,1200]
molarMass = 2* 1.00794# H2 [g/mole]
v_max = 12000
v = np.arange(0,v_max,1)
y = np.array(0)
i = 0
fig,ax = plt.subplots()
while (i<len(T)):
    # maybe store each y, but here they are appended to compare their heights
    y = np.append(y, maxwellBoltzmannWithVelocitySquared(fig,ax,v,molarMass,T[i],crayonForMolecularH2))
    i+=1
# ====================================================================================================
# Include Atomic H
molarMass = 1.00794# H [g/mole]
i=0
while (i<len(T)):
    # maybe store each y, but here they are appended to compare their heights
    y = np.append(y, maxwellBoltzmannWithVelocitySquared(fig,ax,v,molarMass,T[i],crayonForAtomicH))
    i+=1
adjustPlotHeight(ax,max(y))
# ====================================================================================================


# ====================================================================================================
# Escape Velocities
# ====================================================================================================
v_escape = escapeVelocity()
v_escapeMars = escapeVelocityOfMars()
plotEscapeVelocity(fig,ax,v_escape,max(y),'b--',"Escape Velocity of Earth")
plotEscapeVelocity(fig,ax,v_escapeMars,max(y),'r--',"Escape Velocity of Mars")
plotEscapeVelocity(fig3,ax3,v_escape,max(y3),'b--',"Escape Velocity of Earth")
# ====================================================================================================


# ====================================================================================================
plt.show()
