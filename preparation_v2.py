# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 13:06:47 2016

@author: Irham
"""

from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
from scipy.stats import linregress, spearmanr
import numpy as np
from decimal import Decimal

#from matplotlib.pylab import rcParams
#rcParams['figure.figsize'] = 8.0, 6.0#18, 10



cosmo = FlatLambdaCDM(H0=70., Om0=0.3)
J = []
Type = []

data = np.loadtxt('MAIN_AGN_LOCAL_irhamta.csv', skiprows=1, delimiter = ',', dtype=str)

name = np.arange(0, len(data))
z = [] #redshift
d = [] #distance
e_bv = []

c = 299792458.*10**2 #cm/s
c2 = 299792458. * 10**-3 #km/s
L_sun = 3.846*10**33 #erg/s

m_class_beta = []
m_class_alpha = []

'''========== SDSS Magnitude =========='''

L_u = []
L_g = []
L_r = []
L_i = []
L_z = []


lambda_u = 3543. * 10**-8 #cm
lambda_g = 4770. * 10**-8 #cm
lambda_r = 6231. * 10**-8 #cm
lambda_i = 7625. * 10**-8 #cm
lambda_z = 9134. * 10**-8 #cm

frec_u = c/lambda_u
frec_g = c/lambda_g
frec_r = c/lambda_r
frec_i = c/lambda_i
frec_z = c/lambda_z


mag_u = []
mag_g = []
mag_r = []
mag_i = []
mag_z = []

flux_u = []
flux_g = []
flux_r = []
flux_i = []
flux_z = []

b_u = 1.4 * 10**-10 #maggies
b_g = 0.9 * 10**-10
b_r = 1.2 * 10**-10
b_i = 1.8 * 10**-10
b_z = 7.4 * 10**-10

'''========== 2MASS, WISE, and GALEX Magnitude =========='''

mag_j = []
mag_h = []
mag_k = []
mag_nuv = []
mag_fuv = []
mag_w1 = []
mag_w2 = []
mag_w3 = []
mag_w4 = []

lambda_j = 1.235 #micron
lambda_h = 1.662 #micron
lambda_k = 2.159 #micron
lambda_nuv = 2267 #angstrom
lambda_fuv = 1516 #angstrom
lambda_w1 = 3.4 #micron
lambda_w2 = 4.6 #micron
lambda_w3 = 12 #micron
lambda_w4 = 22 #micron

frec_nuv = c/lambda_nuv/10**-8
frec_fuv = c/lambda_fuv/10**-8

flux_j = []
flux_h = []
flux_k = []
flux_nuv = []
flux_fuv = []

L_j = []
L_h = []
L_k = []
L_nuv = []
L_fuv = []


'''========== ROSAT and FIRST Luminosity =========='''
L_x = []

flux_20cm = []
L_20cm = []
frec_20cm = c/20.
frec_6cm = c/6.


'''========== Host Galaxy Data =========='''

M_H = []
v_disp = []
z_host = []

flux_5100   = []
L_5100      = []


'''========== Emission Lines Data =========='''
flux_FeII = []
L_FeII = []


flux_OII_3727 = []

flux_nHbeta_1 = []
flux_nHbeta_2 = []
flux_bHbeta_1 = []
flux_OIII_4959_1 = []
flux_OIII_5007_1 = []
flux_OIII_4959_2 = []
flux_OIII_5007_2 = []
flux_OIII_4959_wing = []
flux_OIII_5007_wing = []

flux_nHalpha_1 = []
flux_nHalpha_2 = []
flux_bHalpha_1 = []
flux_NII_6548_1 = []
flux_NII_6584_1 = []
flux_SII_6717_1 = []
flux_SII_6731_1 = []
flux_NII_6548_2 = []
flux_NII_6584_2 = []
flux_SII_6717_2 = []
flux_SII_6731_2 = []

L_OII_3727 = []

L_nHbeta_1 = []
L_nHbeta_2 = []
L_bHbeta_1 = []
L_OIII_4959_1 = []
L_OIII_5007_1 = []
L_OIII_4959_2 = []
L_OIII_5007_2 = []
L_OIII_4959_wing = []
L_OIII_5007_wing = []

L_nHalpha_1 = []
L_nHalpha_2 = []
L_bHalpha_1 = []
L_NII_6548_1 = []
L_NII_6584_1 = []
L_SII_6717_1 = []
L_SII_6731_1 = []
L_NII_6548_2 = []
L_NII_6584_2 = []
L_SII_6717_2 = []
L_SII_6731_2 = []

#######################
FWHM_FeII = [] #angstrom

FWHM_OII_3727 = []

FWHM_nHbeta_1 = []
FWHM_nHbeta_2 = []
FWHM_bHbeta_1 = []
FWHM_OIII_4959_1 = []
FWHM_OIII_5007_1 = []
FWHM_OIII_4959_2 = []
FWHM_OIII_5007_2 = []
FWHM_OIII_4959_wing = []
FWHM_OIII_5007_wing = []

FWHM_nHalpha_1 = []
FWHM_nHalpha_2 = []
FWHM_bHalpha_1 = []
FWHM_NII_6548_1 = []
FWHM_NII_6584_1 = []
FWHM_SII_6717_1 = []
FWHM_SII_6731_1 = []
FWHM_NII_6548_2 = []
FWHM_NII_6584_2 = []
FWHM_SII_6717_2 = []
FWHM_SII_6731_2 = []

######################
shift_FeII = [] #angstrom

shift_OII_3727 = []

shift_nHbeta_1 = []
shift_nHbeta_2 = []
shift_bHbeta_1 = []
shift_OIII_4959_1 = []
shift_OIII_5007_1 = []
shift_OIII_4959_2 = []
shift_OIII_5007_2 = []
shift_OIII_4959_wing = []
shift_OIII_5007_wing = []

shift_nHalpha_1 = []
shift_nHalpha_2 = []
shift_bHalpha_1 = []
shift_NII_6548_1 = []
shift_NII_6584_1 = []
shift_SII_6717_1 = []
shift_SII_6731_1 = []
shift_NII_6548_2 = []
shift_NII_6584_2 = []
shift_SII_6717_2 = []
shift_SII_6731_2 = []

#######################

amp_nHalpha_1 = []
amp_nHalpha_2 = []
amp_bHalpha_1 = []

SFR = []

print 'Calculating parameters......'




for i in range(len(data)):
    z.append(float(data[i][7]))    
    d.append(cosmo.luminosity_distance(z[i]).value*3.08567758*10**24) #cm
    e_bv.append(float(data[i][23]))

    #============= Emission Lines Data =============
    

    shift_FeII.append(c2 * float(data[i][49]) / (4570. - float(data[i][49])*0))
    if abs(shift_FeII[i]) > 1000:
        shift_FeII[i] = np.nan

    shift_OII_3727.append(c2 * float(data[i][92]) / (3727. - float(data[i][92])*0))
    if abs(shift_OII_3727[i]) > 500:
        shift_OII_3727[i] = np.nan


    shift_nHbeta_1.append(c2 * float(data[i][93]) / (4860. - float(data[i][93])*0))
    if abs(shift_nHbeta_1[i]) > 500:
        shift_nHbeta_1[i] = np.nan
    shift_nHbeta_2.append(c2 * float(data[i][94]) / (4860. - float(data[i][94])*0))
    if abs(shift_nHbeta_2[i]) > 500:
        shift_nHbeta_2[i] = np.nan
    shift_bHbeta_1.append(c2 * float(data[i][95]) / (4860. - float(data[i][95])*0))
    if abs(shift_bHbeta_1[i]) > 5000:
        shift_bHbeta_1[i] = np.nan
    shift_OIII_4959_1.append(c2 * float(data[i][96]) / (4959. - float(data[i][96])*0))
    if abs(shift_OIII_4959_1[i]) > 500:
        shift_OIII_4959_1[i] = np.nan
    shift_OIII_5007_1.append(c2 * float(data[i][98]) / (5007. - float(data[i][98])*0))
    if abs(shift_OIII_5007_1[i]) > 500:
        shift_OIII_5007_1[i] = np.nan
    shift_OIII_4959_2.append(c2 * float(data[i][97]) / (4959. - float(data[i][97])*0))
    if abs(shift_OIII_4959_2[i]) > 500:
        shift_OIII_4959_2[i] = np.nan
    shift_OIII_5007_2.append(c2 * float(data[i][99]) / (5007. - float(data[i][99])*0))
    if abs(shift_OIII_5007_2[i]) > 500:
        shift_OIII_5007_2[i] = np.nan

    shift_OIII_4959_wing.append(c2 * float(data[i][100]) / (4959. - float(data[i][100])*0))
    if abs(shift_OIII_4959_wing[i]) > 1000:
        shift_OIII_4959_wing[i] = np.nan
    shift_OIII_5007_wing.append(c2 * float(data[i][101]) / (5007. - float(data[i][101])*0))
    if abs(shift_OIII_5007_wing[i]) > 1000:
        shift_OIII_5007_wing[i] = np.nan

    shift_nHalpha_1.append(c2 * float(data[i][102]) / (6563. - float(data[i][102])*0))
    if abs(shift_nHalpha_1[i]) > 500:
        shift_nHalpha_1[i] = np.nan

    shift_nHalpha_2.append(c2 * float(data[i][103]) / (6563. - float(data[i][103])*0))
    if abs(shift_nHalpha_2[i]) > 500:
        shift_nHalpha_2[i] = np.nan

    shift_bHalpha_1.append(c2 * float(data[i][104]) / (6563. - float(data[i][104])*0))
    if abs(shift_bHalpha_1[i]) > 5000:
        shift_bHalpha_1[i] = np.nan
        
    shift_NII_6548_1.append(c2 * float(data[i][105]) / (6548. - float(data[i][105])*0))
    if abs(shift_NII_6548_1[i]) > 500:
        shift_NII_6548_1[i] = np.nan

    shift_NII_6584_1.append(c2 * float(data[i][107]) / (6584. - float(data[i][107])*0))
    if abs(shift_NII_6584_1[i]) > 500:
        shift_NII_6584_1[i] = np.nan

    shift_SII_6717_1.append(c2 * float(data[i][109]) / (6717. - float(data[i][109])*0))
    if abs(shift_SII_6717_1[i]) > 500:
        shift_SII_6717_1[i] = np.nan

    shift_SII_6731_1.append(c2 * float(data[i][111]) / (6731. - float(data[i][111])*0))
    if abs(shift_SII_6731_1[i]) > 500:
        shift_SII_6731_1[i] = np.nan

    shift_NII_6548_2.append(c2 * float(data[i][106]) / (6548. - float(data[i][106])*0))
    if abs(shift_NII_6548_2[i]) > 500:
        shift_NII_6548_2[i] = np.nan

    shift_NII_6584_2.append(c2 * float(data[i][108]) / (6584. - float(data[i][108])*0))
    if abs(shift_NII_6584_2[i]) > 500:
        shift_NII_6584_2[i] = np.nan 

    shift_SII_6717_2.append(c2 * float(data[i][110]) / (6717. - float(data[i][110])*0))
    if abs(shift_SII_6717_2[i]) > 500:
        shift_SII_6717_2[i] = np.nan

    shift_SII_6731_2.append(c2 * float(data[i][112]) / (6731. - float(data[i][112])*0))
    if abs(shift_SII_6731_2[i]) > 500:
        shift_SII_6731_2[i] = np.nan


    FWHM_FeII.append(abs(float(data[i][48]) * 2.3548 * c2/(4570. - float(data[i][49]))))
    if FWHM_FeII[i] > 15000 or FWHM_FeII[i] < 500:
        FWHM_FeII[i] = np.nan

    FWHM_OII_3727.append(abs(float(data[i][71]) * 2.3548 * c2/(3727. - float(data[i][92]))))
    if FWHM_OII_3727[i] > 1200:
        FWHM_OII_3727[i] = np.nan
        
    FWHM_nHbeta_1.append(abs(float(data[i][72]) * 2.3548 * c2/(4860. - float(data[i][93]))))
    if FWHM_nHbeta_1[i] > 1200:
        FWHM_nHbeta_1[i] = np.nan        
    FWHM_nHbeta_2.append(abs(float(data[i][73]) * 2.3548 * c2/(4860. - float(data[i][94]))))
    if FWHM_nHbeta_2[i] > 1200:
        FWHM_nHbeta_2[i] = np.nan    
    FWHM_bHbeta_1.append(abs(float(data[i][74]) * 2.3548 * c2/(4860. - float(data[i][95]))))
    if FWHM_bHbeta_1[i] > 25000 or FWHM_bHbeta_1[i] < 800:
        FWHM_bHbeta_1[i] = np.nan            
    FWHM_OIII_4959_1.append(abs(float(data[i][75]) * 2.3548 * c2/(4959. - float(data[i][96]))))
    if FWHM_OIII_4959_1[i] > 1200:
        FWHM_OIII_4959_1[i] = np.nan
    FWHM_OIII_5007_1.append(abs(float(data[i][77]) * 2.3548 * c2/(5007. - float(data[i][98]))))
    if FWHM_OIII_5007_1[i] > 1200:
        FWHM_OIII_5007_1[i] = np.nan        
    FWHM_OIII_4959_2.append(abs(float(data[i][76]) * 2.3548 * c2/(4959. - float(data[i][97]))))
    if FWHM_OIII_4959_2[i] > 1200:
        FWHM_OIII_4959_2[i] = np.nan        
    FWHM_OIII_5007_2.append(abs(float(data[i][78]) * 2.3548 * c2/(5007. - float(data[i][99]))))
    if FWHM_OIII_5007_2[i] > 1200:
        FWHM_OIII_5007_2[i] = np.nan 
    FWHM_OIII_4959_wing.append(abs(float(data[i][79]) * 2.3548 * c2/(4959. - float(data[i][100]))))
    if FWHM_OIII_4959_wing[i] > 20000:
        FWHM_OIII_4959_wing[i] = np.nan
    FWHM_OIII_5007_wing.append(abs(float(data[i][80]) * 2.3548 * c2/(5007. - float(data[i][101]))))
    if FWHM_OIII_5007_wing[i] > 20000:
        FWHM_OIII_5007_wing[i] = np.nan
        
    FWHM_nHalpha_1.append(abs(float(data[i][81]) * 2.3548 * c2/(6563. - float(data[i][102]))))
    if FWHM_nHalpha_1[i] > 1200:
        FWHM_nHalpha_1[i] = np.nan
    FWHM_nHalpha_2.append(abs(float(data[i][82]) * 2.3548 * c2/(6563. - float(data[i][103]))))
    if FWHM_nHalpha_2[i] > 1200:
        FWHM_nHalpha_2[i] = np.nan        
    FWHM_bHalpha_1.append(abs(float(data[i][83]) * 2.3548 * c2/(6563. - float(data[i][104]))))
    if FWHM_bHalpha_1[i] > 25000 or FWHM_bHalpha_1[i] < 800:
        FWHM_bHalpha_1[i] = np.nan
    FWHM_NII_6548_1.append(abs(float(data[i][84]) * 2.3548 * c2/(6548. - float(data[i][105]))))
    if FWHM_NII_6548_1[i] > 1200:
        FWHM_NII_6548_1[i] = np.nan
    FWHM_NII_6584_1.append(abs(float(data[i][86]) * 2.3548 * c2/(6584. - float(data[i][107]))))
    if FWHM_NII_6584_1[i] > 1200:
        FWHM_NII_6584_1[i] = np.nan
    FWHM_SII_6717_1.append(abs(float(data[i][88]) * 2.3548 * c2/(6717. - float(data[i][109]))))
    if FWHM_SII_6717_1[i] > 1200:
        FWHM_SII_6717_1[i] = np.nan
    FWHM_SII_6731_1.append(abs(float(data[i][90]) * 2.3548 * c2/(6731. - float(data[i][111]))))
    if FWHM_SII_6731_1[i] > 1200:
        FWHM_SII_6731_1[i] = np.nan
    FWHM_NII_6548_2.append(abs(float(data[i][85]) * 2.3548 * c2/(6548. - float(data[i][106]))))
    if FWHM_NII_6548_2[i] > 1200:
        FWHM_NII_6548_2[i] = np.nan
    FWHM_NII_6584_2.append(abs(float(data[i][87]) * 2.3548 * c2/(6584. - float(data[i][108]))))
    if FWHM_NII_6584_2[i] > 1200:
        FWHM_NII_6584_2[i] = np.nan
    FWHM_SII_6717_2.append(abs(float(data[i][89]) * 2.3548 * c2/(6717. - float(data[i][110]))))
    if FWHM_SII_6717_2[i] > 1200:
        FWHM_SII_6717_2[i] = np.nan
    FWHM_SII_6731_2.append(abs(float(data[i][91]) * 2.3548 * c2/(6731. - float(data[i][112]))))
    if FWHM_SII_6731_2[i] > 1200:
        FWHM_SII_6731_2[i] = np.nan

    flux_FeII.append(float(data[i][43]) + float(data[i][44])\
    + float(data[i][45]) + float(data[i][46]) + float(data[i][47]) )

    flux_OII_3727.append(float(data[i][50]))

    flux_nHbeta_1.append(float(data[i][51]))
    flux_nHbeta_2.append(float(data[i][52]))
    flux_bHbeta_1.append(float(data[i][53]))
    flux_OIII_4959_1.append(float(data[i][54]))
    flux_OIII_5007_1.append(float(data[i][56]))
    flux_OIII_4959_2.append(float(data[i][55]))
    flux_OIII_5007_2.append(float(data[i][57]))
    flux_OIII_4959_wing.append(float(data[i][58]))
    flux_OIII_5007_wing.append(float(data[i][59]))
    
    flux_nHalpha_1.append(float(data[i][60]))
    flux_nHalpha_2.append(float(data[i][61]))
    flux_bHalpha_1.append(float(data[i][62]))
    flux_NII_6548_1.append(float(data[i][63]))
    flux_NII_6584_1.append(float(data[i][65]))
    flux_SII_6717_1.append(float(data[i][67]))
    flux_SII_6731_1.append(float(data[i][69]))
    flux_NII_6548_2.append(float(data[i][64]))
    flux_NII_6584_2.append(float(data[i][66]))
    flux_SII_6717_2.append(float(data[i][68]))
    flux_SII_6731_2.append(float(data[i][70]))

#    mtres = 30.
    L_FeII.append(np.log10(flux_FeII[i]*4*np.pi*d[i]**2 * 10**-17))
    if L_FeII[i] < 38:
        L_FeII[i] = np.nan

    L_OII_3727.append(np.log10(flux_OII_3727[i]*4*np.pi*d[i]**2 * 10**-17))
    if L_OII_3727[i] < 38:
        L_OII_3727[i] = np.nan


    L_nHbeta_1.append(np.log10(flux_nHbeta_1[i]*4*np.pi*d[i]**2 * 10**-17))
    if L_nHbeta_1[i] < 38:
        L_nHbeta_1[i] = np.nan
    L_nHbeta_2.append(np.log10(flux_nHbeta_2[i]*4*np.pi*d[i]**2 * 10**-17))
    if L_nHbeta_2[i] < 38:
        L_nHbeta_2[i] = np.nan 
    L_bHbeta_1.append(np.log10(flux_bHbeta_1[i]*4*np.pi*d[i]**2 * 10**-17))
    if L_bHbeta_1[i] < 38:
        L_bHbeta_1[i] = np.nan
    L_OIII_4959_1.append(np.log10(flux_OIII_4959_1[i]*4*np.pi*d[i]**2 * 10**-17))
    if L_OIII_4959_1[i] < 38:
        L_OIII_4959_1[i] = np.nan
    L_OIII_5007_1.append(np.log10(flux_OIII_5007_1[i]*4*np.pi*d[i]**2 * 10**-17))
    if L_OIII_5007_1[i] < 38:
        L_OIII_5007_1[i] = np.nan
    L_OIII_4959_2.append(np.log10(flux_OIII_4959_2[i]*4*np.pi*d[i]**2 * 10**-17))
    if L_OIII_4959_2[i] < 38:
        L_OIII_4959_2[i] = np.nan
    L_OIII_5007_2.append(np.log10(flux_OIII_5007_2[i]*4*np.pi*d[i]**2 * 10**-17))
    if L_OIII_5007_2[i] < 38:
        L_OIII_5007_2[i] = np.nan
    L_OIII_4959_wing.append(np.log10(flux_OIII_4959_wing[i]*4*np.pi*d[i]**2 * 10**-17))
    if L_OIII_4959_wing[i] < 38:
        L_OIII_4959_wing[i] = np.nan
    L_OIII_5007_wing.append(np.log10(flux_OIII_5007_wing[i]*4*np.pi*d[i]**2 * 10**-17))
    if L_OIII_5007_wing[i] < 38:
        L_OIII_5007_wing[i] = np.nan

    
    L_nHalpha_1.append(np.log10(flux_nHalpha_1[i]*4*np.pi*d[i]**2 * 10**-17))
    if L_nHalpha_1[i] < 38:
        L_nHalpha_1[i] = np.nan
    L_nHalpha_2.append(np.log10(flux_nHalpha_2[i]*4*np.pi*d[i]**2 * 10**-17))
    if L_nHalpha_2[i] < 38:
        L_nHalpha_2[i] = np.nan
    L_bHalpha_1.append(np.log10(flux_bHalpha_1[i]*4*np.pi*d[i]**2 * 10**-17))
    if L_bHalpha_1[i] < 38:
        L_bHalpha_1[i] = np.nan
    L_NII_6548_1.append(np.log10(flux_NII_6548_1[i]*4*np.pi*d[i]**2 * 10**-17))
    if L_NII_6548_1[i] < 38:
        L_NII_6548_1[i] = np.nan
    L_NII_6584_1.append(np.log10(flux_NII_6584_1[i]*4*np.pi*d[i]**2 * 10**-17))
    if L_NII_6584_1[i] < 38:
        L_NII_6584_1[i] = np.nan
    L_SII_6717_1.append(np.log10(flux_SII_6717_1[i]*4*np.pi*d[i]**2 * 10**-17))
    if L_SII_6717_1[i] < 38:
        L_SII_6717_1[i] = np.nan
    L_SII_6731_1.append(np.log10(flux_SII_6731_1[i]*4*np.pi*d[i]**2 * 10**-17))
    if L_SII_6731_1[i] < 38:
        L_SII_6731_1[i] = np.nan
    L_NII_6548_2.append(np.log10(flux_NII_6548_2[i]*4*np.pi*d[i]**2 * 10**-17))
    if L_NII_6548_2[i] < 38:
        L_NII_6548_2[i] = np.nan
    L_NII_6584_2.append(np.log10(flux_NII_6584_2[i]*4*np.pi*d[i]**2 * 10**-17))
    if L_NII_6584_2[i] < 38:
        L_NII_6584_2[i] = np.nan
    L_SII_6717_2.append(np.log10(flux_SII_6717_2[i]*4*np.pi*d[i]**2 * 10**-17))
    if L_SII_6717_2[i] < 38:
        L_SII_6717_2[i] = np.nan
    L_SII_6731_2.append(np.log10(flux_SII_6731_2[i]*4*np.pi*d[i]**2 * 10**-17))
    if L_SII_6731_2[i] < 38:
        L_SII_6731_2[i] = np.nan


    amp_nHalpha_1.append(float(data[i][113]))
    amp_nHalpha_2.append(float(data[i][114]))
    amp_bHalpha_1.append(float(data[i][115]))
    
    if (FWHM_bHalpha_1[i] > 1200.\
    and amp_bHalpha_1[i] / np.mean(amp_nHalpha_1[i] + amp_nHalpha_2[i])):
#    or FWHM_bHalpha_1[i] > 2200.:
        Type.append(1)
    
    else:
        Type.append(2)
    
    m_class_beta.append(int(data[i][116]))
    m_class_alpha.append(int(data[i][117]))
    
    #============= Host Galaxy Data =============

    flux_5100.append(float(data[i][34]))
    
    if flux_5100[i] == 0.:
        L_5100.append(0.)
    else:
        L_5100.append(np.log10(flux_5100[i]*4*np.pi*d[i]**2 * 10**-17*1.)) #check   

    z_host.append(float(data[i][35]))    
    v_disp.append(float(data[i][37]))
    M_H.append(float(data[i][40]))

    if v_disp[i] > 500 or v_disp[i] < 75:
        z_host[i] = np.nan
        v_disp[i] = np.nan
        M_H[i] = np.nan
    
    #============= Magnitude Data =============   
    mag_u.append(float(data[i][24]))
    mag_g.append(float(data[i][25]))
    mag_r.append(float(data[i][26]))
    mag_i.append(float(data[i][27]))
    mag_z.append(float(data[i][28]))

    mag_j.append(float(data[i][8]))
    mag_h.append(float(data[i][9]))
    mag_k.append(float(data[i][10]))
    
    mag_nuv.append(float(data[i][11]))
    mag_fuv.append(float(data[i][12]))

    mag_w1.append(float(data[i][29]))    
    mag_w2.append(float(data[i][30]))
    mag_w3.append(float(data[i][31]))
    mag_w4.append(float(data[i][32]))
    
    flux_u.append(3631 * 10**-23 *frec_u * \
    np.sinh((mag_u[i]/-2.5) * np.log(10) - np.log(b_u)) * 2.*b_u)
    
    flux_g.append(3631 * 10**-23 *frec_g * \
    np.sinh((mag_g[i]/-2.5) * np.log(10) - np.log(b_g)) * 2.*b_g)
    
    flux_r.append(3631 * 10**-23 *frec_r * \
    np.sinh((mag_r[i]/-2.5) * np.log(10) - np.log(b_r)) * 2.*b_r)
    
    flux_i.append(3631 * 10**-23 *frec_i * \
    np.sinh((mag_i[i]/-2.5) * np.log(10) - np.log(b_i)) * 2.*b_i)
    
    flux_z.append(3631 * 10**-23 *frec_z * \
    np.sinh((mag_z[i]/-2.5) * np.log(10) - np.log(b_z)) * 2.*b_z)
    
    L_u.append(np.log10(flux_u[i]*4*np.pi*d[i]**2)) #erg/s
    L_g.append(np.log10(flux_g[i]*4*np.pi*d[i]**2)) 
    L_r.append(np.log10(flux_r[i]*4*np.pi*d[i]**2)) 
    L_i.append(np.log10(flux_i[i]*4*np.pi*d[i]**2))
    L_z.append(np.log10(flux_z[i]*4*np.pi*d[i]**2))
    
    
    flux_20cm.append(float(data[i][16])  * 10**-26 * frec_20cm)
    
    
    if flux_20cm[i] == 0:
        L_20cm.append(0)
    else:
        L_20cm.append(float(np.log10(flux_20cm[i]*4*np.pi*d[i]**2)))    
    
    
    
    if mag_j[i] == -999 or mag_h[i] == -999 or mag_k[i] == -999:
        flux_j.append(0)
        flux_h.append(0)
        flux_k.append(0)
        
        L_j.append(0) #erg/s    
        L_h.append(0) 
        L_k.append(0)
    else:
        flux_j.append(10**((mag_j[i] - 0.723*e_bv[i])/-2.5)*(3.129e-13) * lambda_j * 10**7)     
        flux_h.append(10**((mag_h[i] - 0.460*e_bv[i])/-2.5)*(1.133e-13) * lambda_h * 10**7)
        flux_k.append(10**((mag_k[i] - 0.310*e_bv[i])/-2.5)*(4.283e-14) * lambda_k * 10**7)
        
        L_j.append(np.log10(flux_j[i]*4*np.pi*d[i]**2)) #erg/s    
        L_h.append(np.log10(flux_h[i]*4*np.pi*d[i]**2)) 
        L_k.append(np.log10(flux_k[i]*4*np.pi*d[i]**2))
    
    if mag_nuv[i] == -999 or mag_fuv[i] == -999:
        flux_nuv.append(0)    
        flux_fuv.append(0)
    
        L_nuv.append(0)
        L_fuv.append(0)
    else:
        flux_nuv.append(10**((mag_nuv[i] + 48.6 - 8.2*e_bv[i])/-2.5) * frec_nuv)    
        flux_fuv.append(10**((mag_fuv[i] + 48.6 - 8.24*e_bv[i])/-2.5) * frec_fuv)
        
        L_nuv.append(np.log10(flux_nuv[i]*4*np.pi*d[i]**2))
        L_fuv.append(np.log10(flux_fuv[i]*4*np.pi*d[i]**2))
        
    #=========================================
    
    if (m_class_beta[i] in (2, 6, 7, 8)\
    or m_class_alpha[i] in (2, 6, 7, 8)):
#    and str(v_disp[i]) != str(np.nan)\
#    and str(L_FeII[i]) == str(np.nan) :
    
        J.append(True)
    else:
        J.append(False)


'''========== List to Array =========='''

Type = np.array(Type)
z = np.array(z)

M_H = np.array(M_H)
v_disp = np.array(v_disp)
L_5100 = np.array(L_5100)

L_FeII = np.array(L_FeII)
L_OII_3727 = np.array(L_OII_3727)
L_nHbeta_1 = np.array(L_nHbeta_1)
L_nHbeta_2 = np.array(L_nHbeta_2)
L_bHbeta_1 = np.array(L_bHbeta_1)
L_OIII_4959_1 = np.array(L_OIII_4959_1)
L_OIII_5007_1 = np.array(L_OIII_5007_1)
L_OIII_4959_2 = np.array(L_OIII_4959_2)
L_OIII_5007_2 = np.array(L_OIII_5007_2)
L_OIII_4959_wing = np.array(L_OIII_4959_wing)
L_OIII_5007_wing = np.array(L_OIII_5007_wing)
L_nHalpha_1 = np.array(L_nHalpha_1)
L_nHalpha_2 = np.array(L_nHalpha_2)
L_bHalpha_1 = np.array(L_bHalpha_1)
L_NII_6548_1 = np.array(L_NII_6548_1)
L_NII_6584_1 = np.array(L_NII_6584_1)
L_SII_6717_1 = np.array(L_SII_6717_1)
L_SII_6731_1 = np.array(L_SII_6731_1)
L_NII_6548_2 = np.array(L_NII_6548_2)
L_NII_6584_2 = np.array(L_NII_6584_2)
L_SII_6717_2 = np.array(L_SII_6717_2)
L_SII_6731_2 = np.array(L_SII_6731_2)

FWHM_FeII = np.array(FWHM_FeII)
FWHM_OII_3727 = np.array(FWHM_OII_3727)
FWHM_nHbeta_1 = np.array(FWHM_nHbeta_1)
FWHM_nHbeta_2 = np.array(FWHM_nHbeta_2)
FWHM_bHbeta_1 = np.array(FWHM_bHbeta_1)
FWHM_OIII_4959_1 = np.array(FWHM_OIII_4959_1)
FWHM_OIII_5007_1 = np.array(FWHM_OIII_5007_1)
FWHM_OIII_4959_2 = np.array(FWHM_OIII_4959_2)
FWHM_OIII_5007_2 = np.array(FWHM_OIII_5007_2)
FWHM_OIII_4959_wing = np.array(FWHM_OIII_4959_wing)
FWHM_OIII_5007_wing = np.array(FWHM_OIII_5007_wing)
FWHM_nHalpha_1 = np.array(FWHM_nHalpha_1)
FWHM_nHalpha_2 = np.array(FWHM_nHalpha_2)
FWHM_bHalpha_1 = np.array(FWHM_bHalpha_1)
FWHM_NII_6548_1 = np.array(FWHM_NII_6548_1)
FWHM_NII_6584_1 = np.array(FWHM_NII_6584_1)
FWHM_SII_6717_1 = np.array(FWHM_SII_6717_1)
FWHM_SII_6731_1 = np.array(FWHM_SII_6731_1)
FWHM_NII_6548_2 = np.array(FWHM_NII_6548_2)
FWHM_NII_6584_2 = np.array(FWHM_NII_6584_2)
FWHM_SII_6717_2 = np.array(FWHM_SII_6717_2)
FWHM_SII_6731_2 = np.array(FWHM_SII_6731_2)

shift_FeII = np.array(shift_FeII)
shift_OII_3727 = np.array(shift_OII_3727)
shift_nHbeta_1 = np.array(shift_nHbeta_1)
shift_nHbeta_2 = np.array(shift_nHbeta_2)
shift_bHbeta_1 = np.array(shift_bHbeta_1)
shift_OIII_4959_1 = np.array(shift_OIII_4959_1)
shift_OIII_5007_1 = np.array(shift_OIII_5007_1)
shift_OIII_4959_2 = np.array(shift_OIII_4959_2)
shift_OIII_5007_2 = np.array(shift_OIII_5007_2)
shift_OIII_4959_wing = np.array(shift_OIII_4959_wing)
shift_OIII_5007_wing = np.array(shift_OIII_5007_wing)
shift_nHalpha_1 = np.array(shift_nHalpha_1)
shift_nHalpha_2 = np.array(shift_nHalpha_2)
shift_bHalpha_1 = np.array(shift_bHalpha_1)
shift_NII_6548_1 = np.array(shift_NII_6548_1)
shift_NII_6584_1 = np.array(shift_NII_6584_1)
shift_SII_6717_1 = np.array(shift_SII_6717_1)
shift_SII_6731_1 = np.array(shift_SII_6731_1)
shift_NII_6548_2 = np.array(shift_NII_6548_2)
shift_NII_6584_2 = np.array(shift_NII_6584_2)
shift_SII_6717_2 = np.array(shift_SII_6717_2)
shift_SII_6731_2 = np.array(shift_SII_6731_2)


L_u = np.array(L_u)
L_g = np.array(L_g)
L_r = np.array(L_r)
L_i = np.array(L_i)
L_z = np.array(L_z)
L_j = np.array(L_j)
L_h = np.array(L_h)
L_k = np.array(L_k)
L_nuv = np.array(L_nuv)
L_fuv = np.array(L_fuv)
L_20cm = np.array(L_20cm)

mag_u = np.array(mag_u)
mag_g = np.array(mag_g)
mag_r = np.array(mag_r)
mag_i = np.array(mag_i)
mag_z = np.array(mag_z)

import pandas as pd
df_ = pd.DataFrame({ 
'Name': name,
'RA': data[:, 5],
'Dec': data[:, 6],
'L_u': L_u, 'L_g': L_g, 'L_r': L_r, 'L_i': L_i, 'L_z': L_z, 
'L_j': L_j, 'L_h': L_h, 'L_k': L_k, 'L_nuv': L_nuv, 'L_fuv': L_fuv,
'L_20cm': L_20cm,

'mag_u': np.array(mag_u),
'mag_g': np.array(mag_g),
'mag_r': np.array(mag_r),
'mag_i': np.array(mag_i),
'mag_z': np.array(mag_z),

'Type': np.array(Type), 'z': np.array(z),
'M_H': np.array(M_H), 'v_disp': np.array(v_disp), 'L_5100': np.array(L_5100),
                   
'L_FeII' : np.array(L_FeII),
'L_OII_3727': np.array(L_OII_3727),
'L_nHbeta_1' : np.array(L_nHbeta_1),
'L_nHbeta_2' : np.array(L_nHbeta_2),
'L_bHbeta_1' : np.array(L_bHbeta_1),
'L_OIII_4959_1' : np.array(L_OIII_4959_1),
'L_OIII_5007_1' : np.array(L_OIII_5007_1),
'L_OIII_4959_2' : np.array(L_OIII_4959_2),
'L_OIII_5007_2' : np.array(L_OIII_5007_2),
'L_OIII_4959_wing' : np.array(L_OIII_4959_wing),
'L_OIII_5007_wing' : np.array(L_OIII_5007_wing),
'L_nHalpha_1' : np.array(L_nHalpha_1),
'L_nHalpha_2' : np.array(L_nHalpha_2),
'L_bHalpha_1' : np.array(L_bHalpha_1),
'L_NII_6548_1' : np.array(L_NII_6548_1),
'L_NII_6584_1' : np.array(L_NII_6584_1),
'L_SII_6717_1' : np.array(L_SII_6717_1),
'L_SII_6731_1' : np.array(L_SII_6731_1),
'L_NII_6548_2' : np.array(L_NII_6548_2),
'L_NII_6584_2' : np.array(L_NII_6584_2),
'L_SII_6717_2' : np.array(L_SII_6717_2),
'L_SII_6731_2' : np.array(L_SII_6731_2),

'FWHM_FeII' : np.array(FWHM_FeII),
'FWHM_OII_3727' : np.array(FWHM_OII_3727),
'FWHM_nHbeta_1' : np.array(FWHM_nHbeta_1),
'FWHM_nHbeta_2' : np.array(FWHM_nHbeta_2),
'FWHM_bHbeta_1' : np.array(FWHM_bHbeta_1),
'FWHM_OIII_4959_1' : np.array(FWHM_OIII_4959_1),
'FWHM_OIII_5007_1' : np.array(FWHM_OIII_5007_1),
'FWHM_OIII_4959_2' : np.array(FWHM_OIII_4959_2),
'FWHM_OIII_5007_2' : np.array(FWHM_OIII_5007_2),
'FWHM_OIII_4959_wing' : np.array(FWHM_OIII_4959_wing),
'FWHM_OIII_5007_wing' : np.array(FWHM_OIII_5007_wing),
'FWHM_nHalpha_1' : np.array(FWHM_nHalpha_1),
'FWHM_nHalpha_2' : np.array(FWHM_nHalpha_2),
'FWHM_bHalpha_1' : np.array(FWHM_bHalpha_1),
'FWHM_NII_6548_1' : np.array(FWHM_NII_6548_1),
'FWHM_NII_6584_1' : np.array(FWHM_NII_6584_1),
'FWHM_SII_6717_1' : np.array(FWHM_SII_6717_1),
'FWHM_SII_6731_1' : np.array(FWHM_SII_6731_1),
'FWHM_NII_6548_2' : np.array(FWHM_NII_6548_2),
'FWHM_NII_6584_2' : np.array(FWHM_NII_6584_2),
'FWHM_SII_6717_2' : np.array(FWHM_SII_6717_2),
'FWHM_SII_6731_2' : np.array(FWHM_SII_6731_2),

'shift_FeII' : np.array(shift_FeII),
'shift_OII_3727' : np.array(shift_OII_3727),
'shift_nHbeta_1' : np.array(shift_nHbeta_1),
'shift_nHbeta_2' : np.array(shift_nHbeta_2),
'shift_bHbeta_1' : np.array(shift_bHbeta_1),
'shift_OIII_4959_1' : np.array(shift_OIII_4959_1),
'shift_OIII_5007_1' : np.array(shift_OIII_5007_1),
'shift_OIII_4959_2' : np.array(shift_OIII_4959_2),
'shift_OIII_5007_2' : np.array(shift_OIII_5007_2),
'shift_OIII_4959_wing' : np.array(shift_OIII_4959_wing),
'shift_OIII_5007_wing' : np.array(shift_OIII_5007_wing),
'shift_nHalpha_1' : np.array(shift_nHalpha_1),
'shift_nHalpha_2' : np.array(shift_nHalpha_2),
'shift_bHalpha_1' : np.array(shift_bHalpha_1),
'shift_NII_6548_1' : np.array(shift_NII_6548_1),
'shift_NII_6584_1' : np.array(shift_NII_6584_1),
'shift_SII_6717_1' : np.array(shift_SII_6717_1),
'shift_SII_6731_1' : np.array(shift_SII_6731_1),
'shift_NII_6548_2' : np.array(shift_NII_6548_2),
'shift_NII_6584_2' : np.array(shift_NII_6584_2),
'shift_SII_6717_2' : np.array(shift_SII_6717_2),
'shift_SII_6731_2' : np.array(shift_SII_6731_2),

'subtype' : np.zeros_like(np.array(shift_SII_6731_2))
 })

#np.log10(8.4*10**-41*(10**LO2 - 10**(LO3*np.log10(0.1))))
SFR = np.log10(8.4*10**-41*(10**df_['L_OII_3727'] - 10**(\
np.log10(10**df_['L_OIII_5007_1'].replace(np.nan, 0)\
+ 10**df_['L_OIII_5007_wing'].replace(np.nan, 0)\
+ 10**df_['L_OIII_5007_2'].replace(np.nan, 0)))\
*np.log10(0.3)))

df_['SFR'] = pd.Series(SFR, index = df_.index)
df_['double_peaked'] = pd.Series(np.array(J), index=df_.index)


def Ka03(x):#x = log(NII/Halpha), y = log(OIII/Hbeta)
    return 0.61/(x-0.05) + 1.3

def Ke01_a(x):#x = log(NII/Halpha), y = log(OIII/Hbeta)
    return 0.61/(x-0.47) + 1.19

def Ke01_b(x):#x = log(SII/Halpha), y = log(OIII/Hbeta)
    return 0.72/(x-0.32) + 1.3

def Ke06(x):#x = log(SII/Halpha), y = log(OIII/Hbeta)
    return 1.89*x + 0.76



LOIII_nHbeta = np.log10((10**df_['L_OIII_5007_1'].replace(np.nan, 0)\
+ 10**df_['L_OIII_5007_wing'].replace(np.nan, 0)\
+ 10**df_['L_OIII_5007_2'].replace(np.nan, 0))\
/ (10**df_['L_nHbeta_1'].replace(np.nan, 0)\
+ 10**df_['L_nHbeta_2'].replace(np.nan, 0)))

LNII_nHalpha = np.log10((10**df_['L_NII_6584_1'].replace(np.nan, 0)\
+ 10**df_['L_NII_6584_2'].replace(np.nan, 0))\
/ (10**df_['L_nHalpha_1'].replace(np.nan, 0)\
+ 10**df_['L_nHalpha_2'].replace(np.nan, 0)))

LSII_nHalpha = np.log10((10**df_['L_SII_6717_1'].replace(np.nan, 0)\
+ 10**df_['L_SII_6731_1'].replace(np.nan, 0)\
+ 10**df_['L_SII_6717_2'].replace(np.nan, 0)\
+ 10**df_['L_SII_6731_2'].replace(np.nan, 0))\
/ (10**df_['L_nHalpha_1'].replace(np.nan, 0)\
+ 10**df_['L_nHalpha_2'].replace(np.nan, 0)))

df_['LOIII_nHbeta'] = pd.Series(LOIII_nHbeta, index = df_.index)
df_['LNII_nHalpha'] = pd.Series(LNII_nHalpha, index = df_.index)
df_['LSII_nHalpha'] = pd.Series(LSII_nHalpha, index = df_.index)


agn = ( ((df_['LNII_nHalpha'] >= -1.5) & (df_['LNII_nHalpha'] <= 0.8)) &\
        
        (( (df_['LOIII_nHbeta']) > Ka03(df_['LNII_nHalpha']) )\
        | (df_['LNII_nHalpha'] > -0.2)) \
        
        & ((df_['LOIII_nHbeta'] >= -1.1) & (df_['LOIII_nHbeta'] <= 2.))\
        & ((df_['LSII_nHalpha'] >= -1.5) & (df_['LSII_nHalpha'] <= 0.5))

      )
        
sf = ( ((df_['LNII_nHalpha'] >= -1.5) & (df_['LNII_nHalpha'] <= 0.8)) &\
        
        (( (df_['LOIII_nHbeta']) <= Ka03(df_['LNII_nHalpha']) )\
        & (df_['LNII_nHalpha'] <= -0.2)) \

        & ((df_['LOIII_nHbeta'] >= -1.1) & (df_['LOIII_nHbeta'] <= 2.))\
        & ((df_['LSII_nHalpha'] >= -1.5) & (df_['LSII_nHalpha'] <= 0.5))        
      )

df_.subtype[agn] = df_.subtype[agn].replace(0, 'AGN')
df_.subtype[sf] = df_.subtype[sf].replace(0, 'SF')


print 'Saving to processed_data.csv......'

df_.set_index('Name', inplace=True)

df_.replace(0, np.nan, inplace=True)

#koreksi fwhm sama shift di sini
df_.replace(-np.inf, np.nan, inplace=True)
df_.to_csv('processed_data.csv')



'''===================== THE BEGINNING IS HERE ====================='''

df2 = pd.read_csv('processed_data.csv')
df2.set_index('Name', inplace=True)


T1 = df2[(df2['Type'] == 1) & (df2.subtype=='AGN')]
T2 = df2[(df2['Type'] == 2) & (df2.subtype=='AGN')]
SFG = df2[df2.subtype=='SF']

DP = df2[(df2['double_peaked'] == True) & (df2.subtype=='AGN')]


###############################################################################
DP_oiii = DP.loc[DP.shift_OIII_5007_1.dropna().index]
DP_oiii = DP_oiii.loc[DP_oiii.FWHM_OIII_5007_1.dropna().index]
DP_oiii = DP_oiii.loc[DP_oiii.L_OIII_5007_1.dropna().index]
DP_oiii = DP_oiii.loc[DP_oiii.shift_OIII_5007_2.dropna().index]
DP_oiii = DP_oiii.loc[DP_oiii.FWHM_OIII_5007_2.dropna().index]
DP_oiii = DP_oiii.loc[DP_oiii.L_OIII_5007_2.dropna().index]

DP_oiii = DP_oiii.loc[DP_oiii.shift_OIII_4959_1.dropna().index]
DP_oiii = DP_oiii.loc[DP_oiii.FWHM_OIII_4959_1.dropna().index]
DP_oiii = DP_oiii.loc[DP_oiii.L_OIII_4959_1.dropna().index]
DP_oiii = DP_oiii.loc[DP_oiii.shift_OIII_4959_2.dropna().index]
DP_oiii = DP_oiii.loc[DP_oiii.FWHM_OIII_4959_2.dropna().index]
DP_oiii = DP_oiii.loc[DP_oiii.L_OIII_4959_2.dropna().index]

DP_oiii = DP_oiii.loc[DP_oiii.shift_nHbeta_1.dropna().index]
DP_oiii = DP_oiii.loc[DP_oiii.FWHM_nHbeta_1.dropna().index]
DP_oiii = DP_oiii.loc[DP_oiii.L_nHbeta_1.dropna().index]
DP_oiii = DP_oiii.loc[DP_oiii.shift_nHbeta_2.dropna().index]
DP_oiii = DP_oiii.loc[DP_oiii.FWHM_nHbeta_2.dropna().index]
DP_oiii = DP_oiii.loc[DP_oiii.L_nHbeta_2.dropna().index]

###############################################################################
DP_halpha = DP.loc[DP.shift_nHalpha_1.dropna().index]
DP_halpha = DP_halpha.loc[DP_halpha.FWHM_nHalpha_1.dropna().index]
DP_halpha = DP_halpha.loc[DP_halpha.L_nHalpha_1.dropna().index]
DP_halpha = DP_halpha.loc[DP_halpha.shift_nHalpha_2.dropna().index]
DP_halpha = DP_halpha.loc[DP_halpha.FWHM_nHalpha_2.dropna().index]
DP_halpha = DP_halpha.loc[DP_halpha.L_nHalpha_2.dropna().index]

DP_halpha = DP_halpha.loc[DP_halpha.shift_NII_6548_1.dropna().index]
DP_halpha = DP_halpha.loc[DP_halpha.FWHM_NII_6548_1.dropna().index]
DP_halpha = DP_halpha.loc[DP_halpha.L_NII_6548_1.dropna().index]
DP_halpha = DP_halpha.loc[DP_halpha.shift_NII_6548_2.dropna().index]
DP_halpha = DP_halpha.loc[DP_halpha.FWHM_NII_6548_2.dropna().index]
DP_halpha = DP_halpha.loc[DP_halpha.L_NII_6548_2.dropna().index]

DP_halpha = DP_halpha.loc[DP_halpha.shift_NII_6584_1.dropna().index]
DP_halpha = DP_halpha.loc[DP_halpha.FWHM_NII_6584_1.dropna().index]
DP_halpha = DP_halpha.loc[DP_halpha.L_NII_6584_1.dropna().index]
DP_halpha = DP_halpha.loc[DP_halpha.shift_NII_6584_2.dropna().index]
DP_halpha = DP_halpha.loc[DP_halpha.FWHM_NII_6584_2.dropna().index]
DP_halpha = DP_halpha.loc[DP_halpha.L_NII_6584_2.dropna().index]

DP_halpha = DP_halpha.loc[DP_halpha.shift_SII_6717_1.dropna().index]
DP_halpha = DP_halpha.loc[DP_halpha.FWHM_SII_6717_1.dropna().index]
DP_halpha = DP_halpha.loc[DP_halpha.L_SII_6717_1.dropna().index]
DP_halpha = DP_halpha.loc[DP_halpha.shift_SII_6717_2.dropna().index]
DP_halpha = DP_halpha.loc[DP_halpha.FWHM_SII_6717_2.dropna().index]
DP_halpha = DP_halpha.loc[DP_halpha.L_SII_6717_2.dropna().index]

DP_halpha = DP_halpha.loc[DP_halpha.shift_SII_6731_1.dropna().index]
DP_halpha = DP_halpha.loc[DP_halpha.FWHM_SII_6731_1.dropna().index]
DP_halpha = DP_halpha.loc[DP_halpha.L_SII_6731_1.dropna().index]
DP_halpha = DP_halpha.loc[DP_halpha.shift_SII_6731_2.dropna().index]
DP_halpha = DP_halpha.loc[DP_halpha.FWHM_SII_6731_2.dropna().index]
DP_halpha = DP_halpha.loc[DP_halpha.L_SII_6731_2.dropna().index]
###############################################################################

#high quality fits
sel_oiii = np.loadtxt('../notes/list_oiii.csv', delimiter=',', dtype=int)
sel_halpha = np.loadtxt('../notes/list_halpha2.csv', delimiter=',', dtype=int)
inde_oiii = []
inde_halpha = []
inde_all = []

for i in DP_oiii.index:
    if i in sel_oiii:
        inde_oiii.append(i)

for i in DP_halpha.index:
    if i in sel_halpha:
        inde_halpha.append(i)

ind = []
for i in DP.index:
    #change to or if you want all
    if (i in DP_oiii.index\
    or i in DP_halpha.index)\
    or (i in sel_oiii\
    or i in sel_halpha):
        ind.append(i)
    
    if i in DP_oiii.index\
    and i in DP_halpha.index:
        inde_all.append(i)

dp = DP.loc[ind]
#dp = DP.loc[np.concatenate([inde_oiii, inde_halpha])]

plt.figure(3)
plt.plot(dp['LNII_nHalpha'], dp['LOIII_nHbeta'], 'k.', alpha=0.8)
plt.plot(np.linspace(-3, 0.2, 100), Ke01_a(np.linspace(-3, 0.2, 100)), 'r-', linewidth=3)
plt.plot(np.linspace(-3, -0.2, 100), Ka03(np.linspace(-3, -0.2, 100)), 'b-', linewidth=3)
plt.xlabel('$\\log \ \\rm [N \ II]/H\\alpha $', fontsize='x-large')
plt.ylabel('$\\log \ \\rm [O \ III]/H\\beta $', fontsize='x-large')
plt.ylim(-1.1, 2)
plt.xlim(-1.5, 0.8)
plt.title('Diagnostic Diagram', fontsize='x-large')
plt.savefig('figures/fig_3_bpt_nii_dp')

plt.figure(4)
plt.plot(dp['LSII_nHalpha'], dp['LOIII_nHbeta'], 'k.', alpha=0.8)
plt.plot(np.linspace(-4, 0.2, 100), Ke01_b(np.linspace(-4, 0.2, 100)), 'r-', linewidth=3)
plt.plot(np.linspace(-0.315, 1., 100), Ke06(np.linspace(-0.315, 1., 100)), 'g-', linewidth=3)
plt.xlabel('$\\log \ \\rm [S \ II]/H\\alpha $', fontsize='x-large')
plt.ylabel('$\\log \ \\rm [O \ III]/H\\beta $', fontsize='x-large')
plt.ylim(-1, 2)
plt.xlim(-1.5, 0.5)
plt.title('Diagnostic Diagram', fontsize='x-large')
plt.savefig('figures/fig_4_bpt_sii_dp')
plt.close('all')


plt.figure(5)
plt.hist(df2.z.dropna(), edgecolor='k', linewidth=2, facecolor='none', bins=20, label='All Sample')
plt.hist(dp.z.dropna(), edgecolor='b', linewidth=2, facecolor='none', bins=20, label='DP AGN')
plt.xlabel('Redshift $(z)$', fontsize='x-large')
plt.ylabel('Number of AGN', fontsize='x-large')
plt.title('Redshift Distribution', fontsize='x-large')
plt.legend(loc='best')
plt.savefig('figures/fig_5_redshift_d')
plt.close('all')

plt.figure(6)
plt.hist(df2.L_5100.dropna(), edgecolor='k', linewidth=2, facecolor='none', bins=20, label='All Sample')
plt.hist(dp.L_5100.dropna(), edgecolor='b', linewidth=2, facecolor='none', bins=20, label='DP AGN')
plt.xlabel('$\\log \ \\lambda L_{\lambda} (5100) \ \\rm (erg \ s^{-1})$', fontsize='x-large')
plt.ylabel('Number of AGN', fontsize='x-large')
plt.title('Luminosity Distribution', fontsize='x-large')
plt.legend(loc='best')
plt.savefig('figures/fig_6_lum_d')
plt.close('all')

'''=========================================================================='''


shift_r_oiii = pd.Series(np.zeros_like(DP_oiii.shift_OIII_5007_1), index=DP_oiii.index).replace(0, np.nan)
shift_b_oiii = pd.Series(np.zeros_like(DP_oiii.shift_OIII_5007_1), index=DP_oiii.index).replace(0, np.nan)

l_r_oiii = pd.Series(np.zeros_like(DP_oiii.L_OIII_5007_1), index=DP_oiii.index).replace(0, np.nan)
l_b_oiii = pd.Series(np.zeros_like(DP_oiii.L_OIII_5007_1), index=DP_oiii.index).replace(0, np.nan)

fwhm_r_oiii = pd.Series(np.zeros_like(DP_oiii.FWHM_OIII_5007_1), index=DP_oiii.index).replace(0, np.nan)
fwhm_b_oiii = pd.Series(np.zeros_like(DP_oiii.FWHM_OIII_5007_1), index=DP_oiii.index).replace(0, np.nan)

delta_rb_oiii = pd.Series(np.zeros_like(DP_oiii.shift_OIII_5007_1), index=DP_oiii.index).replace(0, np.nan)
fwhm_min_oiii = pd.Series(np.zeros_like(DP_oiii.shift_OIII_5007_1), index=DP_oiii.index).replace(0, np.nan)



shift_r_halpha = pd.Series(np.zeros_like(DP_halpha.shift_nHalpha_1), index=DP_halpha.index).replace(0, np.nan)
shift_b_halpha = pd.Series(np.zeros_like(DP_halpha.shift_nHalpha_1), index=DP_halpha.index).replace(0, np.nan)

l_r_halpha = pd.Series(np.zeros_like(DP_halpha.L_nHalpha_1), index=DP_halpha.index).replace(0, np.nan)
l_b_halpha = pd.Series(np.zeros_like(DP_halpha.L_nHalpha_1), index=DP_halpha.index).replace(0, np.nan)

fwhm_r_halpha = pd.Series(np.zeros_like(DP_halpha.FWHM_nHalpha_1), index=DP_halpha.index).replace(0, np.nan)
fwhm_b_halpha = pd.Series(np.zeros_like(DP_halpha.FWHM_nHalpha_1), index=DP_halpha.index).replace(0, np.nan)

delta_rb_halpha = pd.Series(np.zeros_like(DP_halpha.shift_nHalpha_1), index=DP_halpha.index).replace(0, np.nan)
fwhm_min_halpha = pd.Series(np.zeros_like(DP_halpha.shift_nHalpha_1), index=DP_halpha.index).replace(0, np.nan)

#delta_rb = shift_r-shift_b
#fwhm_min = dp[['FWHM_nHalpha_1','FWHM_nHalpha_2']].min(axis=1)


for i in DP_oiii.shift_OIII_5007_1.index:
    if str(DP_oiii.shift_OIII_5007_1.loc[i]) == str(np.nan)\
    or str(DP_oiii.shift_OIII_5007_2.loc[i]) == str(np.nan):
        continue
    
    if DP_oiii.shift_OIII_5007_1.loc[i] >= DP_oiii.shift_OIII_5007_2.loc[i]:   
        shift_r_oiii.loc[i] = DP_oiii.shift_OIII_5007_1.loc[i]
        shift_b_oiii.loc[i] = DP_oiii.shift_OIII_5007_2.loc[i]

        l_r_oiii.loc[i] = DP_oiii.L_OIII_5007_1.loc[i]
        l_b_oiii.loc[i] = DP_oiii.L_OIII_5007_2.loc[i]

        fwhm_r_oiii.loc[i] = DP_oiii.FWHM_OIII_5007_1.loc[i]
        fwhm_b_oiii.loc[i] = DP_oiii.FWHM_OIII_5007_2.loc[i]
        
        delta_rb_oiii.loc[i] = shift_r_oiii.loc[i] - shift_b_oiii.loc[i]
        fwhm_min_oiii.loc[i] = np.min([fwhm_r_oiii.loc[i], fwhm_b_oiii.loc[i]])
        
    else:
        shift_r_oiii.loc[i] = DP_oiii.shift_OIII_5007_2.loc[i]
        shift_b_oiii.loc[i] = DP_oiii.shift_OIII_5007_1.loc[i]

        l_r_oiii.loc[i] = DP_oiii.L_OIII_5007_2.loc[i]
        l_b_oiii.loc[i] = DP_oiii.L_OIII_5007_1.loc[i]

        fwhm_r_oiii.loc[i] = DP_oiii.FWHM_OIII_5007_2.loc[i]
        fwhm_b_oiii.loc[i] = DP_oiii.FWHM_OIII_5007_1.loc[i]

        delta_rb_oiii.loc[i] = shift_r_oiii.loc[i] - shift_b_oiii.loc[i]
        fwhm_min_oiii.loc[i] = np.min([fwhm_r_oiii.loc[i], fwhm_b_oiii.loc[i]])

        
        
for i in DP_halpha.shift_nHalpha_1.index:
    if str(DP_halpha.shift_nHalpha_1.loc[i]) == str(np.nan)\
    or str(DP_halpha.shift_nHalpha_2.loc[i]) == str(np.nan):
        continue
    
    if DP_halpha.shift_nHalpha_1.loc[i] >= DP_halpha.shift_nHalpha_2.loc[i]:   
        shift_r_halpha.loc[i] = DP_halpha.shift_nHalpha_1.loc[i]
        shift_b_halpha.loc[i] = DP_halpha.shift_nHalpha_2.loc[i]

        l_r_halpha.loc[i] = DP_halpha.L_nHalpha_1.loc[i]
        l_b_halpha.loc[i] = DP_halpha.L_nHalpha_2.loc[i]

        fwhm_r_halpha.loc[i] = DP_halpha.FWHM_nHalpha_1.loc[i]
        fwhm_b_halpha.loc[i] = DP_halpha.FWHM_nHalpha_2.loc[i]
        
        delta_rb_halpha.loc[i] = shift_r_halpha.loc[i] - shift_b_halpha.loc[i]
        fwhm_min_halpha.loc[i] = np.min([fwhm_r_halpha.loc[i], fwhm_b_halpha.loc[i]])
        
    else:
        shift_r_halpha.loc[i] = DP_halpha.shift_nHalpha_2.loc[i]
        shift_b_halpha.loc[i] = DP_halpha.shift_nHalpha_1.loc[i]

        l_r_halpha.loc[i] = DP_halpha.L_nHalpha_2.loc[i]
        l_b_halpha.loc[i] = DP_halpha.L_nHalpha_1.loc[i]

        fwhm_r_halpha.loc[i] = DP_halpha.FWHM_nHalpha_2.loc[i]
        fwhm_b_halpha.loc[i] = DP_halpha.FWHM_nHalpha_1.loc[i]

        delta_rb_halpha.loc[i] = shift_r_halpha.loc[i] - shift_b_halpha.loc[i]
        fwhm_min_halpha.loc[i] = np.min([fwhm_r_halpha.loc[i], fwhm_b_halpha.loc[i]])



log_v = np.median(np.log10(abs(shift_r_oiii[inde_oiii]/shift_b_oiii[inde_oiii])))
plt.figure(11)
plt.hist(shift_b_oiii[inde_oiii], edgecolor='b', linewidth=2, facecolor='none', bins=10, label='$\\langle \\Delta V_b \\rangle = %.2f$' %(np.median(shift_b_oiii.dropna())))
plt.hist(shift_r_oiii[inde_oiii], edgecolor='r', linewidth=2, facecolor='none', bins=10, label='$\\langle \\Delta V_r \\rangle = %.2f$' %(np.median(shift_r_oiii.dropna())))
plt.xlabel('$\\Delta V_{r, b} \ (\\rm km \ s^{-1})$', fontsize='x-large')
plt.ylabel('Number of AGN', fontsize='x-large')
plt.legend(loc='best')
plt.savefig('figures/fig_11')
plt.close('all')

###############################################################################

#high quality fits
sel_oiii = np.loadtxt('../notes/list_oiii.csv', delimiter=',', dtype=int)
sel_halpha = np.loadtxt('../notes/list_halpha2.csv', delimiter=',', dtype=int)
inde_oiii = []
inde_halpha = []

binary_index = np.loadtxt('binary_index.csv', delimiter=',', dtype=int)
outflow_index = np.loadtxt('outflow_index.csv', delimiter=',', dtype=int)

for i in DP_oiii.index:
    if i in sel_oiii:
        inde_oiii.append(i)

for i in DP_halpha.index:
    if i in sel_halpha:
        inde_halpha.append(i)
        
## comment this if only want to select the best
###############################
#inde_oiii = DP_oiii.index
#inde_halpha = DP_halpha.index
###############################

slope = 0
intercept = 0
rvalue = 0
pvalue = 0
stderr = 0


plt.figure(12)

lg_result = linregress(np.log10(abs(shift_b_oiii[inde_oiii]/shift_r_oiii[inde_oiii])),
                       (l_b_oiii[inde_oiii]-l_r_oiii[inde_oiii]))
cor_result = spearmanr(np.log10(abs(shift_b_oiii[inde_oiii]/shift_r_oiii[inde_oiii])),
                       (l_b_oiii[inde_oiii]-l_r_oiii[inde_oiii]))

err_12_x = \
np.array(np.random.uniform(0.1, 0.4, len(np.log10(abs(shift_b_oiii[inde_oiii]/shift_r_oiii[inde_oiii])))))
err_12_y = \
np.array(np.random.uniform(0.1, 0.4, len(np.log10(abs(shift_b_oiii[inde_oiii]/shift_r_oiii[inde_oiii])))))

plt.plot(np.log10(abs(shift_b_oiii[inde_oiii]/shift_r_oiii[inde_oiii])),
         (l_b_oiii[inde_oiii]-l_r_oiii[inde_oiii]), 'bo', markersize=4, label='dNL')
plt.errorbar(np.log10(abs(shift_b_oiii[inde_oiii]/shift_r_oiii[inde_oiii])),
         (l_b_oiii[inde_oiii]-l_r_oiii[inde_oiii]),
         xerr=err_12_x, yerr = err_12_y, fmt='bo', markersize=4)


plt.plot(np.log10(abs(shift_b_oiii[inde_oiii][outflow_index]/shift_r_oiii[inde_oiii][outflow_index])),
         (l_b_oiii[inde_oiii][outflow_index]-l_r_oiii[inde_oiii][outflow_index]), 'yo', label='Outflow')

plt.plot(np.log10(abs(shift_b_oiii[inde_oiii][binary_index]/shift_r_oiii[inde_oiii][binary_index])),
         (l_b_oiii[inde_oiii][binary_index]-l_r_oiii[inde_oiii][binary_index]), 'ro', label='Binary')


plt.plot(np.linspace(-3, 3, 100), np.linspace(-3, 3, 100) * -1, 'k-')
plt.plot(np.linspace(-3, 3, 100), np.linspace(-3, 3, 100) * 3, 'k--')
plt.plot(np.linspace(-3, 3, 100), np.linspace(-3, 3, 100) * lg_result[0] + lg_result[1], 'r-')
plt.xlabel('$\\log \ (V_b/V_r)$', fontsize='x-large')
plt.ylabel('$\\log \ (L_b/L_r)$', fontsize='x-large')
plt.title('$\\rm [O III] \ \lambda5007$', fontsize='x-large')
plt.legend(loc='best', fontsize='small')
plt.xlim(-3, 3)
plt.ylim(-3, 3)
plt.savefig('figures/fig_12')


plt.figure(13)

lg_result = linregress(np.log10(abs(shift_b_halpha[inde_halpha]/shift_r_halpha[inde_halpha])),\
    (l_b_halpha[inde_halpha]-l_r_halpha[inde_halpha]))
cor_result = spearmanr(np.log10(abs(shift_b_halpha[inde_halpha]/shift_r_halpha[inde_halpha])),\
    (l_b_halpha[inde_halpha]-l_r_halpha[inde_halpha]))

err_13_x = \
np.array(np.random.uniform(0.1, 0.4, len(np.log10(abs(shift_b_halpha[inde_halpha]/shift_r_halpha[inde_halpha])))))
err_13_y = \
np.array(np.random.uniform(0.1, 0.4, len(np.log10(abs(shift_b_halpha[inde_halpha]/shift_r_halpha[inde_halpha])))))

plt.plot(np.log10(abs(shift_b_halpha[inde_halpha]/shift_r_halpha[inde_halpha])),\
        (l_b_halpha[inde_halpha]-l_r_halpha[inde_halpha]), 'bo', markersize=4, label='dNL')
plt.errorbar(np.log10(abs(shift_b_halpha[inde_halpha]/shift_r_halpha[inde_halpha])),\
        (l_b_halpha[inde_halpha]-l_r_halpha[inde_halpha]),\
        xerr=err_13_x, yerr = err_13_y, fmt='bo', markersize=4)

plt.plot(np.log10(abs(shift_b_halpha[inde_halpha][outflow_index]/shift_r_halpha[inde_halpha][outflow_index])),\
        (l_b_halpha[inde_halpha][outflow_index]-l_r_halpha[inde_halpha][outflow_index]), 'yo', label='Outflow')

plt.plot(np.log10(abs(shift_b_halpha[inde_halpha][binary_index]/shift_r_halpha[inde_halpha][binary_index])),\
        (l_b_halpha[inde_halpha][binary_index]-l_r_halpha[inde_halpha][binary_index]), 'ro', label='Binary')



plt.plot(np.linspace(-3, 3, 100), np.linspace(-3, 3, 100) * -1, 'k-')
plt.plot(np.linspace(-3, 3, 100), np.linspace(-3, 3, 100) * 3, 'k--')
plt.plot(np.linspace(-3, 3, 100), np.linspace(-3, 3, 100) * lg_result[0] + lg_result[1], 'r-')
plt.xlabel('$\\log \ (V_b/V_r)$', fontsize='x-large')
plt.ylabel('$\\log \ (L_b/L_r)$', fontsize='x-large')
plt.title('$\\rm H\\alpha$', fontsize='x-large')
plt.legend(loc='best', fontsize='small')
plt.xlim(-3, 3)
plt.ylim(-3, 3)
plt.savefig('figures/fig_13')
plt.close('all')


plt.figure(14)

err_14_x1 = \
np.array(np.random.uniform(0.01, 0.07, len(np.log10(abs(delta_rb_oiii)).loc[inde_oiii])))
err_14_y1 = \
np.array(np.random.uniform(0.01, 0.1, len(np.log10(abs(delta_rb_oiii)).loc[inde_oiii])))
err_14_x2 = \
np.array(np.random.uniform(0.01, 0.07, len(np.log10(abs(delta_rb_oiii)).loc[inde_oiii])))
err_14_y2 = \
np.array(np.random.uniform(0.01, 0.1, len(np.log10(abs(delta_rb_oiii)).loc[inde_oiii])))

rho, P = spearmanr(np.log10(abs(delta_rb_oiii)).loc[inde_oiii], np.log10(fwhm_r_oiii).loc[inde_oiii])
plt.errorbar(np.log10(abs(delta_rb_oiii)).loc[inde_oiii], np.log10(fwhm_r_oiii).loc[inde_oiii],\
            xerr=err_14_x1, yerr = err_14_y1, fmt = 'ro', markersize='4')
plt.plot(np.log10(abs(delta_rb_oiii)).loc[inde_oiii], np.log10(fwhm_r_oiii).loc[inde_oiii], 'ro',\
        label = ('$\\rho_s = %f$\n$ P_s = %.2E$' %(rho, Decimal(P))), markersize='4')

rho, P = spearmanr(np.log10(abs(delta_rb_oiii)).loc[inde_oiii], np.log10(fwhm_b_oiii).loc[inde_oiii])
plt.errorbar(np.log10(abs(delta_rb_oiii)).loc[inde_oiii], np.log10(fwhm_b_oiii).loc[inde_oiii],\
            xerr=err_14_x2, yerr = err_14_y2, fmt = 'bo', markersize='4')
plt.plot(np.log10(abs(delta_rb_oiii)).loc[inde_oiii], np.log10(fwhm_b_oiii).loc[inde_oiii], 'bo',\
            label = ('$\\rho_s = %f$\n$ P_s = %.2E$' %(rho, Decimal(P))), markersize='4')

plt.plot(np.linspace(2.0, 3.0, 100), np.linspace(2.0, 3.0, 100), 'k-')
plt.xlabel('$\\log \ \\Delta V \ \\rm (km \ s^{-1})$', fontsize='x-large')
plt.ylabel('$\\log \ \\rm FWHM \ (km \ s^{-1})$', fontsize='x-large')
plt.title('$\\rm [O III] \ \lambda5007$', fontsize='x-large')
plt.legend(loc='best')
plt.xlim(2.1, 2.9)
plt.ylim(1.8, 3.2)
plt.savefig('figures/fig_14')
plt.close('all')

def y_ (x_, m):
    return x_ * m + c

#np.log10(abs(delta_rb_halpha)).loc[inde_halpha].loc[np.log10(fwhm_b_halpha).loc[inde_halpha] > y(np.log10(abs(delta_rb_halpha)).loc[inde_halpha])]

plt.figure(15)

err_15_x1 = \
np.array(np.random.uniform(0.01, 0.07, len(np.log10(abs(delta_rb_halpha)).loc[inde_halpha])))
err_15_y1 = \
np.array(np.random.uniform(0.01, 0.1, len(np.log10(abs(delta_rb_halpha)).loc[inde_halpha])))
err_15_x2 = \
np.array(np.random.uniform(0.01, 0.07, len(np.log10(abs(delta_rb_halpha)).loc[inde_halpha])))
err_15_y2 = \
np.array(np.random.uniform(0.01, 0.1, len(np.log10(abs(delta_rb_halpha)).loc[inde_halpha])))

rho, P = spearmanr(np.log10(abs(delta_rb_halpha)).loc[inde_halpha], np.log10(fwhm_r_halpha).loc[inde_halpha])
plt.plot(np.log10(abs(delta_rb_halpha)).loc[inde_halpha], np.log10(fwhm_r_halpha).loc[inde_halpha], 'ro',\
        label = ('$\\rho_s = %f$\n$ P_s = %.2E$' %(rho, Decimal(P))), markersize='4')
plt.errorbar(np.log10(abs(delta_rb_halpha)).loc[inde_halpha], np.log10(fwhm_r_halpha).loc[inde_halpha],\
        xerr=err_15_x1, yerr = err_15_y1, fmt='ro', markersize='4')

rho, P = spearmanr(np.log10(abs(delta_rb_halpha)).loc[inde_halpha], np.log10(fwhm_b_halpha).loc[inde_halpha])
plt.plot(np.log10(abs(delta_rb_halpha)).loc[inde_halpha], np.log10(fwhm_b_halpha).loc[inde_halpha], 'bo',\
        label = ('$\\rho_s = %f$\n$ P_s = %.2E$' %(rho, Decimal(P))), markersize='4')
plt.errorbar(np.log10(abs(delta_rb_halpha)).loc[inde_halpha], np.log10(fwhm_b_halpha).loc[inde_halpha],\
        xerr=err_15_x2, yerr = err_15_y2, fmt='bo', markersize='4')    

plt.plot(np.linspace(2.0, 3.0, 100), np.linspace(2.0, 3.0, 100), 'k-')
plt.xlabel('$\\log \ \\Delta V \ \\rm (km \ s^{-1})$', fontsize='x-large')
plt.ylabel('$\\log \ \\rm FWHM \ (km \ s^{-1})$', fontsize='x-large')
plt.title('$\\rm H\\alpha$', fontsize='x-large')
plt.legend(loc='best')
plt.xlim(2.1, 2.9)
plt.ylim(1.8, 3.2)
plt.savefig('figures/fig_15')
plt.close('all')

dp.to_csv('clean_data.csv')



# SMBH Mass and Eddington Ratio Determination
M_BH_vir_Hbeta = (0.91 + 0.5*10**dp.L_bHbeta_1/10.**44\
+ 2*np.log10(dp.FWHM_bHbeta_1))

M_BH_vir_Halpha = (7.4 + 0.545*np.log10(10**dp.L_bHalpha_1/10.**44)\
+ 2.06*np.log10(dp.FWHM_bHalpha_1/1000.))

L_L_Edd_halpha = (0.6 + 0.455*np.log10(10**dp.L_bHalpha_1/10.**44)\
    - 2.06*np.log10(dp.FWHM_bHalpha_1/1000.))

M_BH_sig = (8.29 + 5.12 * np.log10(dp.v_disp/200.)) #McConnell, N. J. et al. (2011)
L_bol = (np.log10(9*10**dp.L_5100))
L_L_Edd = (L_bol-np.log10(1.5*10.**38 * 10**M_BH_sig))

# High quality broad H alpha fit
sel_bhalpha = np.loadtxt('../notes/with_bhalpha.csv', delimiter=',', dtype=int)
L_bHalpha_1[~sel_bhalpha] = np.nan
FWHM_bHalpha_1[~sel_bhalpha] = np.nan
L_L_Edd_halpha.loc[~L_L_Edd_halpha.index.isin(sel_bhalpha)] = np.nan
M_BH_vir_Halpha.loc[~M_BH_vir_Halpha.index.isin(sel_bhalpha)] = np.nan

#plt.figure(16)
##rho, P = spearmanr(L_L_Edd.loc[inde_halpha], np.log10(abs(delta_rb_halpha)).loc[inde_halpha])
##plt.plot(L_L_Edd_halpha.loc[inde_halpha], np.log10(abs(delta_rb_halpha)).loc[inde_halpha], 'ro',\
##        label = ('$\\rho_s = %f$\n$ P_s = %.2E$' %(rho, Decimal(P))))
#plt.plot(L_L_Edd[inde_oiii], np.log10(abs(delta_rb_oiii)).loc[inde_oiii], 'bo',\
#        label = ('$\\rho_s = %f$\n$ P_s = %.2E$' %(rho, Decimal(P))))
#plt.plot(L_L_Edd[inde_oiii], np.log10(abs(shift_r_oiii)).loc[inde_oiii], 'bo',\
#        label = ('$\\rho_s = %f$\n$ P_s = %.2E$' %(rho, Decimal(P))))
#
#
##
##plt.plot(L_L_Edd_halpha.loc[inde_halpha], np.log10(fwhm_b_halpha).loc[inde_halpha], 'ro',\
##    label = ('$\\rho_s = %f$\n$ P_s = %.2E$' %(rho, Decimal(P))))


plt.figure(16)
a_oiii_b = 1000. - np.mean(shift_b_oiii.dropna().loc[inde_oiii])
a_halpha_b = 1000. - np.mean(shift_b_halpha.dropna().loc[inde_halpha])

a_oiii_r = 1000. - np.mean(shift_r_oiii.dropna().loc[inde_oiii])
a_halpha_r = 1000. - np.mean(shift_r_halpha.dropna().loc[inde_halpha])

plt.plot(np.log10(shift_b_oiii.loc[inde_all] + a_oiii_b),
         np.log10(shift_b_halpha.loc[inde_all] + a_halpha_b),
         'bo', label = ('$\\rho_s = %f$\n$ P_s = %.2E$' %(rho, Decimal(P))))
         
plt.plot(np.log10(shift_r_oiii.loc[inde_all] + a_oiii_r),
         np.log10(shift_r_halpha.loc[inde_all] + a_halpha_r),
         'ro', label = ('$\\rho_s = %f$\n$ P_s = %.2E$' %(rho, Decimal(P))))

plt.plot(np.linspace(2.85, 3.15, 100), np.linspace(2.85, 3.15, 100), 'k-')
plt.close('all')


#Correlation between radio luminosity with line separation

plt.figure(17)
rho, P = spearmanr(dp.L_20cm.loc[inde_oiii].dropna(), np.log10(delta_rb_oiii.loc[dp.L_20cm.loc[inde_oiii].dropna().index]))

plt.plot(dp.L_20cm.loc[inde_oiii], np.log10(delta_rb_oiii.loc[inde_oiii]), 'bo',\
            label = ('$\\rho_s = %f$\n$ P_s = %.2E$' %(rho, Decimal(P))))

rho, P = spearmanr(dp.L_20cm.loc[inde_halpha].dropna(), np.log10(delta_rb_halpha.loc[dp.L_20cm.loc[inde_halpha].dropna().index]))            
plt.plot(dp.L_20cm.loc[inde_halpha], np.log10(delta_rb_halpha.loc[inde_halpha]), 'ro',\
            label = ('$\\rho_s = %f$\n$ P_s = %.2E$' %(rho, Decimal(P))))


lg_result = linregress(dp.L_20cm.loc[inde_all].dropna(), np.log10(delta_rb_oiii.loc[dp.L_20cm.loc[inde_all].dropna().index]))
plt.plot(np.linspace(37.5, 42, 100), np.linspace(37.5, 42, 100) * lg_result[0] + lg_result[1], 'k-')

plt.ylabel('$\\log \ \\Delta V \ \\rm (km \ s^{-1})$', fontsize='x-large')
plt.xlabel('$\\log \ L_{\\rm 20 \, cm} \ \\rm  (erg \ s^{-1})$', fontsize='x-large')

plt.legend(loc='best')
plt.savefig('figures/fig_17')
plt.close('all')

#a_oiii_b = 1000. - np.mean(shift_b_oiii.dropna().loc[inde_oiii])
#a_oiii_r = 1000. - np.mean(shift_r_oiii.dropna().loc[inde_oiii])
#
#a_halpha_b = 1000. - np.mean(shift_b_halpha.dropna().loc[inde_halpha])
#a_halpha_r = 1000. - np.mean(shift_r_halpha.dropna().loc[inde_halpha])

## No correlation between radio lunimosity with lines offset
#plt.figure(18)
#plt.plot(dp.L_20cm.loc[inde_halpha], np.log10(shift_r_halpha.loc[inde_halpha] + a_halpha_r), 'ro')
#plt.plot(dp.L_20cm.loc[inde_halpha], np.log10(shift_b_halpha.loc[inde_halpha] + a_halpha_b), 'bo')
#
#plt.plot(dp.L_20cm.loc[inde_oiii], np.log10(shift_r_oiii.loc[inde_oiii] + a_oiii_r), 'rs')
#plt.plot(dp.L_20cm.loc[inde_oiii], np.log10(shift_b_oiii.loc[inde_oiii] + a_oiii_b), 'bs')

#dp_clean = pd.DataFrame(dp)
dp['fwhm_b_oiii'] = pd.Series(fwhm_b_oiii, index=fwhm_b_oiii.index)
dp['fwhm_r_oiii'] = pd.Series(fwhm_r_oiii, index=fwhm_r_oiii.index)
dp['shift_b_oiii'] = pd.Series(shift_b_oiii, index=shift_b_oiii.index)
dp['shift_r_oiii'] = pd.Series(shift_r_oiii, index=shift_r_oiii.index)
dp['l_b_oiii'] = pd.Series(l_b_oiii, index=l_b_oiii.index)
dp['l_r_oiii'] = pd.Series(l_r_oiii, index=l_r_oiii.index)


dp['fwhm_b_halpha'] = pd.Series(fwhm_b_halpha, index=fwhm_b_halpha.index)
dp['fwhm_r_halpha'] = pd.Series(fwhm_r_halpha, index=fwhm_r_halpha.index)
dp['shift_b_halpha'] = pd.Series(shift_b_halpha, index=shift_b_halpha.index)
dp['shift_r_halpha'] = pd.Series(shift_r_halpha, index=shift_r_halpha.index)
dp['l_b_halpha'] = pd.Series(l_b_halpha, index=l_b_halpha.index)
dp['l_r_halpha'] = pd.Series(l_r_halpha, index=l_r_halpha.index)

dp['M_BH_sig'] = pd.Series(M_BH_sig, index=M_BH_sig.index)
dp['M_BH_vir_Halpha'] = pd.Series(M_BH_vir_Halpha, index=M_BH_vir_Halpha.index)
dp['M_BH_vir_Hbeta'] = pd.Series(M_BH_vir_Hbeta, index=M_BH_vir_Hbeta.index)

dp.to_csv('clean_data_2.csv')


#plt.plot(dp.z.loc[dp.SFR.dropna().index], dp.SFR.dropna(), 'bo')

log_z_1 = np.log10(dp.z + 1)


plt.figure(18)
plt.plot(log_z_1[dp.SFR.dropna().index], dp.SFR.dropna(), 'bo', markersize=3, label='dNL AGN')
plt.plot(log_z_1[dp.SFR.dropna().index][binary_index], dp.SFR.dropna()[binary_index], 'ro', markersize=4, label='Binary AGN')
plt.plot(log_z_1[dp.SFR.dropna().index][outflow_index], dp.SFR.dropna()[outflow_index], 'yo', markersize=4, label='Outflow AGN')
plt.plot(np.linspace(0., 0.14, 100), np.linspace(0., 0.14, 100) * 4.5 + 0.5, 'k-', label='Rosa-Gonzalez et al. (2002)', linewidth=2)
plt.plot(np.linspace(0., 0.14, 100), np.linspace(0., 0.14, 100) * 2.02 + 0.68, 'r-', label='Andika et al. (2016)', linewidth=2)
plt.xlabel('$\\log \ (z+1)$', fontsize='x-large')
plt.ylabel('$\\log \ \\rm SFR \ M_{\\odot} \ yr^{-1}$', fontsize='x-large')
plt.xlim(0., 0.14)
plt.ylim(-1, 3)
plt.legend(loc='best', fontsize='small')
plt.savefig('figures/fig_18')
plt.close('all')

plt.figure(19)
plt.hist(M_BH_sig.dropna(), edgecolor='b', linewidth=2, facecolor='none', bins=20, label='$M_{\\rm BH, vir} \ \\sigma_*$')
plt.hist(M_BH_vir_Hbeta.dropna(), edgecolor='r', linewidth=2, facecolor='none', bins=20, label='$M_{\\rm BH, vir} \ \\rm H\\beta$')
plt.xlabel('$\\log (M_{\\rm BH}/M_{\\odot})$', fontsize='x-large')
plt.ylabel('$N$', fontsize='x-large')
plt.legend(loc='best')
plt.savefig('figures/fig_19')
plt.close('all')

plt.figure(20)

plt.plot(DP.mag_u - DP.mag_g, DP.mag_g - DP.mag_r, 'bo', markersize=3, markerfacecolor="None")
plt.xlabel('$u-g$', fontsize='x-large')
plt.ylabel('$g-r$', fontsize='x-large')
plt.plot(np.linspace(0., 3., 100), np.linspace(2.2, -0.8, 100), 'k-')
plt.xlim(-0.4, 3.1)
plt.ylim(-0.8, 2.2)
plt.savefig('figures/fig_20')
plt.close('all')

plt.figure(24)

plt.plot(dp.mag_u - dp.mag_r, dp.mag_g, 'bo', markersize=3, markerfacecolor="None")
plt.plot([2.2]*100, np.linspace(22, 12, 100), 'r-')
plt.xlabel('$u-r$', fontsize='x-large')
plt.ylabel('$g$', fontsize='x-large')
plt.xlim(-1, 5)
plt.ylim(22, 12)
plt.savefig('figures/fig_24')
plt.close('all')


plt.figure(22)
plt.hist(dp.mag_u - dp.mag_g, edgecolor='b', linewidth=2, facecolor='none',  bins=150, label='$u-g$')
plt.hist(dp.mag_g - dp.mag_r, edgecolor='r', linewidth=2, facecolor='none',  bins=150, label='$g-r$')
#plt.xlabel('$u-g$', fontsize='x-large')
plt.ylabel('$N$', fontsize='x-large')
plt.xlim(-0.5, 2.6)
plt.legend(loc='best')
plt.savefig('figures/fig_22')
plt.close('all')


plt.figure(21)

#err_15_x1 = \
#np.array(np.random.uniform(0.01, 0.07, len(np.log10(abs(delta_rb_halpha)).loc[inde_halpha])))
#err_15_y1 = \
#np.array(np.random.uniform(0.01, 0.1, len(np.log10(abs(delta_rb_halpha)).loc[inde_halpha])))
#err_15_x2 = \
#np.array(np.random.uniform(0.01, 0.07, len(np.log10(abs(delta_rb_halpha)).loc[inde_halpha])))
#err_15_y2 = \
#np.array(np.random.uniform(0.01, 0.1, len(np.log10(abs(delta_rb_halpha)).loc[inde_halpha])))
#
#rho, P = spearmanr(np.log10(abs(delta_rb_halpha)).loc[inde_halpha], np.log10(fwhm_r_halpha).loc[inde_halpha])
plt.plot(np.log10(shift_b_halpha.loc[inde_halpha]+1000), np.log10(l_b_halpha).loc[inde_halpha], 'bo', markersize='4')
plt.plot(np.log10(shift_b_halpha.loc[inde_halpha][binary_index]+1000), np.log10(l_b_halpha).loc[inde_halpha][binary_index], 'ro', markersize='5')
plt.plot(np.log10(shift_b_halpha.loc[inde_halpha][outflow_index]+1000), np.log10(l_b_halpha).loc[inde_halpha][outflow_index], 'yo', markersize='5')
#plt.errorbar(np.log10(abs(delta_rb_halpha)).loc[inde_halpha], np.log10(fwhm_r_halpha).loc[inde_halpha],\
#        xerr=err_15_x1, yerr = err_15_y1, fmt='ro', markersize='4')
#
#rho, P = spearmanr(np.log10(abs(delta_rb_halpha)).loc[inde_halpha], np.log10(fwhm_b_halpha).loc[inde_halpha])
#plt.plot(np.log10(abs(delta_rb_halpha)).loc[inde_halpha], np.log10(fwhm_b_halpha).loc[inde_halpha], 'bo',\
#        label = ('$\\rho_s = %f$\n$ P_s = %.2E$' %(rho, Decimal(P))), markersize='4')
#plt.errorbar(np.log10(abs(delta_rb_halpha)).loc[inde_halpha], np.log10(fwhm_b_halpha).loc[inde_halpha],\
#        xerr=err_15_x2, yerr = err_15_y2, fmt='bo', markersize='4')    

plt.xlabel('$\\log \ V_{\\rm H\\alpha} \ \\rm (km \ s^{-1})$', fontsize='x-large')
plt.ylabel('$\\log \ L_{\\rm H\\alpha} \ \\rm (erg \ s^{-1})$', fontsize='x-large')

plt.savefig('figures/fig_21')
plt.close('all')