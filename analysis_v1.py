# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 11:13:52 2016

@author: irham
"""

from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
from scipy.stats import linregress, spearmanr
import numpy as np


cosmo = FlatLambdaCDM(H0=70., Om0=0.3)
J = []
Type = []

data = np.loadtxt('BINARY_AGN.csv', skiprows=2, delimiter = ',', dtype=str)

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



for i in range(len(data)):
    z.append(float(data[i][7]))    
    d.append(cosmo.luminosity_distance(z[i]).value*3.08567758*10**24) #cm
    e_bv.append(float(data[i][23]))

    #============= Emission Lines Data =============

    shift_FeII.append(c2 * float(data[i][49]) / (4570. - float(data[i][49])))

    shift_OII_3727.append(c2 * float(data[i][92]) / (3727. - float(data[i][92])))

    shift_nHbeta_1.append(c2 * float(data[i][93]) / (4860. - float(data[i][93])))
    shift_nHbeta_2.append(c2 * float(data[i][94]) / (4860. - float(data[i][94])))
    shift_bHbeta_1.append(c2 * float(data[i][95]) / (4860. - float(data[i][95])))
    shift_OIII_4959_1.append(c2 * float(data[i][96]) / (4959. - float(data[i][96])))
    shift_OIII_5007_1.append(c2 * float(data[i][98]) / (5007. - float(data[i][98])))
    shift_OIII_4959_2.append(c2 * float(data[i][97]) / (4959. - float(data[i][97])))
    shift_OIII_5007_2.append(c2 * float(data[i][99]) / (5007. - float(data[i][99])))
    shift_OIII_4959_wing.append(c2 * float(data[i][100]) / (4959. - float(data[i][100])))
    shift_OIII_5007_wing.append(c2 * float(data[i][101]) / (5007. - float(data[i][101])))
    
    shift_nHalpha_1.append(c2 * float(data[i][102]) / (6563. - float(data[i][102])))
    shift_nHalpha_2.append(c2 * float(data[i][103]) / (6563. - float(data[i][103])))
    shift_bHalpha_1.append(c2 * float(data[i][104]) / (6563. - float(data[i][104])))
    shift_NII_6548_1.append(c2 * float(data[i][105]) / (6548. - float(data[i][105])))
    shift_NII_6584_1.append(c2 * float(data[i][107]) / (6584. - float(data[i][107])))
    shift_SII_6717_1.append(c2 * float(data[i][109]) / (6717. - float(data[i][109])))
    shift_SII_6731_1.append(c2 * float(data[i][111]) / (6731. - float(data[i][111])))
    shift_NII_6548_2.append(c2 * float(data[i][106]) / (6548. - float(data[i][106])))
    shift_NII_6584_2.append(c2 * float(data[i][108]) / (6584. - float(data[i][108])))
    shift_SII_6717_2.append(c2 * float(data[i][110]) / (6717. - float(data[i][110])))
    shift_SII_6731_2.append(c2 * float(data[i][112]) / (6731. - float(data[i][112])))

    FWHM_FeII.append(float(data[i][48]) * 2.3548 * c2/(4570. - float(data[i][49])))

    FWHM_OII_3727.append(float(data[i][71]) * 2.3548 * c2/(3727. - float(data[i][92])))

    FWHM_nHbeta_1.append(float(data[i][72]) * 2.3548 * c2/(4860. - float(data[i][93])))
    FWHM_nHbeta_2.append(float(data[i][73]) * 2.3548 * c2/(4860. - float(data[i][94])))
    FWHM_bHbeta_1.append(float(data[i][74]) * 2.3548 * c2/(4860. - float(data[i][95])))
    FWHM_OIII_4959_1.append(float(data[i][75]) * 2.3548 * c2/(4959. - float(data[i][96])))
    FWHM_OIII_5007_1.append(float(data[i][77]) * 2.3548 * c2/(5007. - float(data[i][98])))
    FWHM_OIII_4959_2.append(float(data[i][76]) * 2.3548 * c2/(4959. - float(data[i][97])))
    FWHM_OIII_5007_2.append(float(data[i][78]) * 2.3548 * c2/(5007. - float(data[i][99])))
    FWHM_OIII_4959_wing.append(float(data[i][79]) * 2.3548 * c2/(4959. - float(data[i][100])))
    FWHM_OIII_5007_wing.append(float(data[i][80]) * 2.3548 * c2/(5007. - float(data[i][101])))
    
    FWHM_nHalpha_1.append(float(data[i][81]) * 2.3548 * c2/(6563. - float(data[i][102])))
    FWHM_nHalpha_2.append(float(data[i][82]) * 2.3548 * c2/(6563. - float(data[i][103])))
    FWHM_bHalpha_1.append(float(data[i][83]) * 2.3548 * c2/(6563. - float(data[i][104])))
    FWHM_NII_6548_1.append(float(data[i][84]) * 2.3548 * c2/(6548. - float(data[i][105])))
    FWHM_NII_6584_1.append(float(data[i][86]) * 2.3548 * c2/(6584. - float(data[i][107])))
    FWHM_SII_6717_1.append(float(data[i][88]) * 2.3548 * c2/(6717. - float(data[i][109])))
    FWHM_SII_6731_1.append(float(data[i][90]) * 2.3548 * c2/(6731. - float(data[i][111])))
    FWHM_NII_6548_2.append(float(data[i][85]) * 2.3548 * c2/(6548. - float(data[i][106])))
    FWHM_NII_6584_2.append(float(data[i][87]) * 2.3548 * c2/(6584. - float(data[i][108])))
    FWHM_SII_6717_2.append(float(data[i][89]) * 2.3548 * c2/(6717. - float(data[i][110])))
    FWHM_SII_6731_2.append(float(data[i][91]) * 2.3548 * c2/(6731. - float(data[i][112])))





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

    L_FeII.append(np.log10(flux_FeII[i]*4*np.pi*d[i]**2 * 10**-17))

    L_OII_3727.append(np.log10(flux_OII_3727[i]*4*np.pi*d[i]**2 * 10**-17))

    L_nHbeta_1.append(np.log10(flux_nHbeta_1[i]*4*np.pi*d[i]**2 * 10**-17))
    L_nHbeta_2.append(np.log10(flux_nHbeta_2[i]*4*np.pi*d[i]**2 * 10**-17))
    L_bHbeta_1.append(np.log10(flux_bHbeta_1[i]*4*np.pi*d[i]**2 * 10**-17))
    L_OIII_4959_1.append(np.log10(flux_OIII_4959_1[i]*4*np.pi*d[i]**2 * 10**-17))
    L_OIII_5007_1.append(np.log10(flux_OIII_5007_1[i]*4*np.pi*d[i]**2 * 10**-17))
    L_OIII_4959_2.append(np.log10(flux_OIII_4959_2[i]*4*np.pi*d[i]**2 * 10**-17))
    L_OIII_5007_2.append(np.log10(flux_OIII_5007_2[i]*4*np.pi*d[i]**2 * 10**-17))
    L_OIII_4959_wing.append(np.log10(flux_OIII_4959_wing[i]*4*np.pi*d[i]**2 * 10**-17))
    L_OIII_5007_wing.append(np.log10(flux_OIII_5007_wing[i]*4*np.pi*d[i]**2 * 10**-17))
    
    L_nHalpha_1.append(np.log10(flux_nHalpha_1[i]*4*np.pi*d[i]**2 * 10**-17))
    L_nHalpha_2.append(np.log10(flux_nHalpha_2[i]*4*np.pi*d[i]**2 * 10**-17))
    L_bHalpha_1.append(np.log10(flux_bHalpha_1[i]*4*np.pi*d[i]**2 * 10**-17))
    L_NII_6548_1.append(np.log10(flux_NII_6548_1[i]*4*np.pi*d[i]**2 * 10**-17))
    L_NII_6584_1.append(np.log10(flux_NII_6584_1[i]*4*np.pi*d[i]**2 * 10**-17))
    L_SII_6717_1.append(np.log10(flux_SII_6717_1[i]*4*np.pi*d[i]**2 * 10**-17))
    L_SII_6731_1.append(np.log10(flux_SII_6731_1[i]*4*np.pi*d[i]**2 * 10**-17))
    L_NII_6548_2.append(np.log10(flux_NII_6548_2[i]*4*np.pi*d[i]**2 * 10**-17))
    L_NII_6584_2.append(np.log10(flux_NII_6584_2[i]*4*np.pi*d[i]**2 * 10**-17))
    L_SII_6717_2.append(np.log10(flux_SII_6717_2[i]*4*np.pi*d[i]**2 * 10**-17))
    L_SII_6731_2.append(np.log10(flux_SII_6731_2[i]*4*np.pi*d[i]**2 * 10**-17))



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
    or m_class_alpha[i] in (2, 6, 7, 8))\
    and v_disp[i] < 500.\
    and L_FeII[i] < 1. :
    
        J.append(i)


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