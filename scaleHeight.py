import numpy as np
import matplotlib.pyplot as plt

#
#        ____________
#       /           /
#      /           /|
#     /___________/_|_________ H of ligher atom
#     |          |  |
#     |          |  |
#     |          |  |
#     |          |  |
#     |          |  |
#     |          |  |
#     |          |  |             _.-'''''-_
#     |          |  |          _.'          '__-''-_ 
#     |          |  |       .-'      cloud         _\
#     |          |  |       '._                  /'
#     |          |  |          '._-\   _-'''\__-'
#     |          |  |               '-'     \ 
#     |          |  |                \ \  \  \
#     |          |  |              \    \  \  \
#     |          |  |                  \  \   \  rain
#     |          |  |
#     |          |  |                         
#     |          |  |
#     |  ________|__|
#     | /        |  |
#     |/         | /|
#     |__________|/_|__________ H of heavier molecule
#     |          |  |
#     |          |  |
#     |          |  |
#     |          |  |
#     |          |  |
#     |          |  |
#     |          |  |
#     |          |  |
#     |  ________|__|
#     | /        |  /
#     |/         | /
#     |__________|/ .__________ Person, not an ant
#
def H(m,T):
    # parameter: m  - Average atomic mass
    # parameter: T  - Temperature in Kelvin
    NAvogadro = 6.02 * 10**(23) # (Halliday Resnick and Krane, 5e) 
    k = 1.38 * 10**(-23)        # (Halliday Resnick and Krane, 5e)
    m_SI = m * 10**(-3)         # kg
    g = 9.8                     # G*Me*1kg / Re^2 (might get a slightly different 3rd significant figure, but 9.8 to two significant figures is consisent)
    numerator = k*T
    denominator = m_SI / NAvogadro * g
    return(numerator/denominator)

def plotH(ax,H,molecule,crayon):
    ax.fill([0,0,10,10],[0,H,H,0],alpha=0.3,c=crayon)
    ax.text(0,H,molecule)
    return()

mO2 = 32
mH2 = 2
T = 288
HO2 = H(mO2,T)
HO  = H(mO2/2,T)
HH2 = H(mH2,T)

fig,ax = plt.subplots()
plotH(ax,HO2/10**3,"O2",'c')
plotH(ax,HO/10**3,"O",'b')
plotH(ax,HH2/10**3,"H2",'#aabbcc')

title  = "Scale Heights ( H = kT/mg )\nfor Selected Molecules and Atoms\nin a Gravitationally Separated Atmosphere"
ylabel = "Height (km)"
xlabel = "None (Arbitrary)"
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_title(title)
fig.set_figheight(10)
fig.set_figwidth(5)


plt.show()
