# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 22:34:32 2015

@author: irham
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.modeling.models import custom_model
from scipy import interpolate, integrate


o = open('calculation/Fe_II_flux.csv', 'w')
o.write('flux_F\tflux_G\tflux_P\tflux_S\tflux_Z\tsigma_Fe\tshift_Fe\n')
o.close()


#J = np.loadtxt('notes/skipped_j_opt_fix.csv', skiprows=1, dtype=int)
f = open('calculation/continuum_flux_1.csv', 'w')
f.write('F_5100\n')
f.close()

q = 0
'''========== Fe II Template and Spectra Data =========='''

F = np.loadtxt('FeII_template_4000_5500/FeII_rel_int/F_rel_int.txt')
G = np.loadtxt('FeII_template_4000_5500/FeII_rel_int/G_rel_int.txt')
P = np.loadtxt('FeII_template_4000_5500/FeII_rel_int/P_rel_int.txt')
S = np.loadtxt('FeII_template_4000_5500/FeII_rel_int/S_rel_int.txt')
Z = np.loadtxt('FeII_template_4000_5500/FeII_rel_int/IZw1_rel_int.txt')
#agn = np.loadtxt('shifted_redcor/out_spec-0542-51993-0560.fits')

flux_F = []
flux_G = []
flux_P = []
flux_S = []
flux_Z = []
sigma_fe = []
shift_fe = []

wl_cont = []
fl_cont = []

sample = np.loadtxt('MAIN_AGN_LOCAL_irhamta.csv', skiprows=1, dtype=str, delimiter=',')

plate = []
mjd = []
fiberid = []
z = []
e_bv = []

skip = np.loadtxt('notes/skip.csv', delimiter=',', skiprows=1)

for j in range(len(sample)):
    if len(sample[j][2]) < 4:
        plate.append(str(0)+sample[j][2])
    else:
        plate.append(sample[j][2])
    
    mjd.append(sample[j][3])

    if len(sample[j][4]) == 4:
        fiberid.append(sample[j][4])    
    elif len(sample[j][4]) == 3:
        fiberid.append(str(0) + sample[j][4])
    elif len(sample[j][4]) == 2:
        fiberid.append(str(0) + str(0) + sample[j][4])
    else:
        fiberid.append(str(0) + str(0) + str(0) + sample[j][4])
    
    if j < 23781:
        print 'skipped, j = ', j
        continue    
    
    
    if j in skip[:, 1] or j == skip[0, 4] or j not in skip[:, 5]: #or j < 3455:#for until redshift 3.5
        f = open('calculation/continuum_flux_1.csv', 'a')
        #for i in range(len(fl_cont)):
        f.write(str(0) + '\n')

        f.close()

        o = open('calculation/Fe_II_flux.csv', 'a')
        o.write(str(0) + '\t' + str(0) + '\t' + str(0)\
        + '\t' + str(0) + '\t' + str(0) + '\t' + str(0)\
        + '\t' + str(0)+ '\n')
    
        o.close()
        
        print 'skipped, j = ', j
        continue
    
    agn = np.loadtxt('shifted_redcor/out_spec-' + plate[j] + '-' + mjd[j] + '-' + fiberid[j] + '.fits')
    


#'''========== Fitting Near Hbeta ==========='''

    wl_hbeta_cont = []
    fl_hbeta_cont = []
    wl_hbeta = []
    fl_hbeta = []


    wl_halpha_cont = []
    fl_halpha_cont = []
    wl_halpha = []
    fl_halpha = []


    wl_OII_cont = []
    fl_OII_cont = []
    wl_OII = []
    fl_OII = []

    for i in range(len(agn)):#Hbeta window
        if agn[i][0] >= 5095 and agn[i][0] <=5105:
            wl_cont.append(agn[i][0])
            
    
    
        if agn[i][0] >= 4435 and agn[i][0] <=5535:
            wl_hbeta.append(agn[i][0])
            fl_hbeta.append(agn[i][1])
        elif agn[i][0] >= 6000 and agn[i][0] <=7000:
            wl_halpha.append(agn[i][0])
            fl_halpha.append(agn[i][1])
        elif agn[i][0] >= 3650 and agn[i][0] <=3810:#3747:
            wl_OII.append(agn[i][0])
            fl_OII.append(agn[i][1])
        else:
            continue

    #if min(wl_OII) > 3717:
    #    J.append(j)
    #    print 'skipped, j = ', j
    #    continue


    for i in range(len(agn)):#Hbeta continuum window 
#        if agn[i][0] >= 4435 and agn[i][0] <=4700\
#        or agn[i][0] >= 5100 and agn[i][0] <=5535:
        if agn[i][0] >= 4150 and agn[i][0] <= 4250\
        or agn[i][0] >= 4450 and agn[i][0] <= 4600\
        or agn[i][0] >= 5100 and agn[i][0] <= 5500:
            wl_hbeta_cont.append(agn[i][0])
            fl_hbeta_cont.append(agn[i][1])
 

        
        elif agn[i][0] >= 6000 and agn[i][0] <=6250\
        or agn[i][0] >= 6800 and agn[i][0] <=7000: 
            wl_halpha_cont.append(agn[i][0])
            fl_halpha_cont.append(agn[i][1])
    
        elif agn[i][0] >= 3650 and agn[i][0] <= 3717\
        or agn[i][0] >= 3790 and agn[i][0] <= 3810:
            wl_OII_cont.append(agn[i][0])
            fl_OII_cont.append(agn[i][1])
        
        else:
            continue

    @custom_model   
    def power_law_hbeta(x, amplitude = max(fl_hbeta_cont), alpha = 1.5):
        return amplitude * (x**-alpha)

    @custom_model
    def FeII_template_hbeta(x, mF = np.mean(fl_hbeta_cont), \
    mG = np.mean(fl_hbeta_cont), mP = np.mean(fl_hbeta_cont), \
    mS = np.mean(fl_hbeta_cont), mZ = np.mean(fl_hbeta_cont), \
    shift = 0, sigma = 1):
    
        return (abs(mF) * \
        (F[0][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[0][0] + shift) / sigma)**2.)\
        + F[1][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[1][0] + shift) / sigma)**2.)\
        + F[2][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[2][0] + shift) / sigma)**2.)\
        + F[3][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[3][0] + shift) / sigma)**2.)\
        + F[4][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[4][0] + shift) / sigma)**2.)\
        + F[5][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[5][0] + shift) / sigma)**2.)\
        + F[6][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[6][0] + shift) / sigma)**2.)\
        + F[7][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[7][0] + shift) / sigma)**2.)\
        + F[8][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[8][0] + shift) / sigma)**2.)\
        + F[9][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[9][0] + shift) / sigma)**2.)\
        + F[10][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[10][0] + shift) / sigma)**2.)\
        + F[11][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[11][0] + shift) / sigma)**2.)\
        + F[12][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[12][0] + shift) / sigma)**2.)\
        + F[13][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[13][0] + shift) / sigma)**2.)\
        + F[14][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[14][0] + shift) / sigma)**2.)\
        + F[15][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[15][0] + shift) / sigma)**2.)\
        + F[16][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[16][0] + shift) / sigma)**2.)\
        + F[17][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[17][0] + shift) / sigma)**2.)\
        + F[18][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[18][0] + shift) / sigma)**2.)\
        + F[19][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[19][0] + shift) / sigma)**2.)\
        + F[20][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[20][0] + shift) / sigma)**2.))\
    
        + abs(mG) * \
        (+ G[0][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - G[0][0] + shift) / sigma)**2.)\
        + G[1][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - G[1][0] + shift) / sigma)**2.)\
        + G[2][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - G[2][0] + shift) / sigma)**2.)\
        + G[3][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - G[3][0] + shift) / sigma)**2.)\
        + G[4][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - G[4][0] + shift) / sigma)**2.)\
        + G[5][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - G[5][0] + shift) / sigma)**2.)\
        + G[6][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - G[6][0] + shift) / sigma)**2.)\
        + G[7][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - G[7][0] + shift) / sigma)**2.)\
        + G[8][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - G[8][0] + shift) / sigma)**2.)\
        + G[9][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - G[9][0] + shift) / sigma)**2.)\
        + G[10][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - G[10][0] + shift) / sigma)**2.))\
    
        + abs(mP) * \
        (+ P[0][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[0][0] + shift) / sigma)**2.)\
        + P[1][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[1][0] + shift) / sigma)**2.)\
        + P[2][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[2][0] + shift) / sigma)**2.)\
        + P[3][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[3][0] + shift) / sigma)**2.)\
        + P[4][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[4][0] + shift) / sigma)**2.)\
        + P[5][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[5][0] + shift) / sigma)**2.)\
        + P[6][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[6][0] + shift) / sigma)**2.)\
        + P[7][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[7][0] + shift) / sigma)**2.)\
        + P[8][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[8][0] + shift) / sigma)**2.)\
        + P[9][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[9][0] + shift) / sigma)**2.)\
        + P[10][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[10][0] + shift) / sigma)**2.)\
        + P[11][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[11][0] + shift) / sigma)**2.)\
        + P[12][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[12][0] + shift) / sigma)**2.)\
        + P[13][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[13][0] + shift) / sigma)**2.)\
        + P[14][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[14][0] + shift) / sigma)**2.))\
        
        + abs(mS) * \
        (+ S[0][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - S[0][0] + shift) / sigma)**2.)\
        + S[1][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - S[1][0] + shift) / sigma)**2.)\
        + S[2][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - S[2][0] + shift) / sigma)**2.)\
        + S[3][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - S[3][0] + shift) / sigma)**2.)\
        + S[4][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - S[4][0] + shift) / sigma)**2.)\
        + S[5][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - S[5][0] + shift) / sigma)**2.)\
        + S[6][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - S[6][0] + shift) / sigma)**2.))\
    
        + abs(mZ) * \
        (+ Z[0][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[0][0] + shift) / sigma)**2.)\
        + Z[1][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[1][0] + shift) / sigma)**2.)\
        + Z[2][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[2][0] + shift) / sigma)**2.)\
        + Z[3][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[3][0] + shift) / sigma)**2.)\
        + Z[4][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[4][0] + shift) / sigma)**2.)\
        + Z[5][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[5][0] + shift) / sigma)**2.)\
        + Z[6][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[6][0] + shift) / sigma)**2.)\
        + Z[7][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[7][0] + shift) / sigma)**2.)\
        + Z[8][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[8][0] + shift) / sigma)**2.)\
        + Z[9][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[9][0] + shift) / sigma)**2.)\
        + Z[10][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[10][0] + shift) / sigma)**2.)\
        + Z[11][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[11][0] + shift) / sigma)**2.)\
        + Z[12][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[12][0] + shift) / sigma)**2.)\
        + Z[13][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[13][0] + shift) / sigma)**2.)\
        + Z[14][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[14][0] + shift) / sigma)**2.)))

    m_init = FeII_template_hbeta() + power_law_hbeta() #For Hbeta
    fit = fitting.LevMarLSQFitter()
    m = fit(m_init, wl_hbeta_cont, fl_hbeta_cont, maxiter=10**3)


    m_init2 = FeII_template_hbeta(mF = np.mean(fl_halpha_cont), mG = np.mean(fl_halpha_cont), \
    mP = np.mean(fl_halpha_cont), mS = np.mean(fl_halpha_cont), \
    mZ = np.mean(fl_halpha_cont), shift = 0, sigma = 1) \
    + power_law_hbeta(amplitude = min(fl_halpha_cont), alpha = 1.5)
 
    fit = fitting.LevMarLSQFitter()
    m2 = fit(m_init2, wl_halpha_cont, fl_halpha_cont, maxiter=10**3)

#'''
#m_init3 = FeII_template_hbeta(mF = np.mean(fl_OII_cont), mG = np.mean(fl_OII_cont), \
#mP = np.mean(fl_OII_cont), mS = np.mean(fl_OII_cont), \
#mZ = np.mean(fl_OII_cont), shift = 0, sigma = 1) \
#+ power_law_hbeta(amplitude = min(fl_OII_cont), alpha = 1.5)

#fit = fitting.LevMarLSQFitter()
#m3 = fit(m_init3, wl_OII_cont, fl_OII_cont, maxiter=300000)
#'''



#'''========== Fitting Decompotition =========='''
    def F_fun(x, mF, shift, sigma):
        return (abs(mF) * \
        (F[0][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[0][0] + shift) / sigma)**2.)\
        + F[1][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[1][0] + shift) / sigma)**2.)\
        + F[2][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[2][0] + shift) / sigma)**2.)\
        + F[3][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[3][0] + shift) / sigma)**2.)\
        + F[4][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[4][0] + shift) / sigma)**2.)\
        + F[5][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[5][0] + shift) / sigma)**2.)\
        + F[6][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[6][0] + shift) / sigma)**2.)\
        + F[7][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[7][0] + shift) / sigma)**2.)\
        + F[8][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[8][0] + shift) / sigma)**2.)\
        + F[9][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[9][0] + shift) / sigma)**2.)\
        + F[10][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[10][0] + shift) / sigma)**2.)\
        + F[11][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[11][0] + shift) / sigma)**2.)\
        + F[12][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[12][0] + shift) / sigma)**2.)\
        + F[13][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[13][0] + shift) / sigma)**2.)\
        + F[14][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[14][0] + shift) / sigma)**2.)\
        + F[15][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[15][0] + shift) / sigma)**2.)\
        + F[16][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[16][0] + shift) / sigma)**2.)\
        + F[17][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[17][0] + shift) / sigma)**2.)\
        + F[18][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[18][0] + shift) / sigma)**2.)\
        + F[19][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[19][0] + shift) / sigma)**2.)\
        + F[20][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - F[20][0] + shift) / sigma)**2.)))

    def G_fun(x, mG, shift, sigma):
        return (abs(mG) * \
        (+ G[0][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - G[0][0] + shift) / sigma)**2.)\
        + G[1][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - G[1][0] + shift) / sigma)**2.)\
        + G[2][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - G[2][0] + shift) / sigma)**2.)\
        + G[3][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - G[3][0] + shift) / sigma)**2.)\
        + G[4][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - G[4][0] + shift) / sigma)**2.)\
        + G[5][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - G[5][0] + shift) / sigma)**2.)\
        + G[6][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - G[6][0] + shift) / sigma)**2.)\
        + G[7][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - G[7][0] + shift) / sigma)**2.)\
        + G[8][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - G[8][0] + shift) / sigma)**2.)\
        + G[9][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - G[9][0] + shift) / sigma)**2.)\
        + G[10][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - G[10][0] + shift) / sigma)**2.)))


    def P_fun(x, mP, shift, sigma):
        return (abs(mP) * \
        (+ P[0][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[0][0] + shift) / sigma)**2.)\
        + P[1][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[1][0] + shift) / sigma)**2.)\
        + P[2][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[2][0] + shift) / sigma)**2.)\
        + P[3][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[3][0] + shift) / sigma)**2.)\
        + P[4][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[4][0] + shift) / sigma)**2.)\
        + P[5][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[5][0] + shift) / sigma)**2.)\
        + P[6][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[6][0] + shift) / sigma)**2.)\
        + P[7][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[7][0] + shift) / sigma)**2.)\
        + P[8][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[8][0] + shift) / sigma)**2.)\
        + P[9][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[9][0] + shift) / sigma)**2.)\
        + P[10][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[10][0] + shift) / sigma)**2.)\
        + P[11][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[11][0] + shift) / sigma)**2.)\
        + P[12][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[12][0] + shift) / sigma)**2.)\
        + P[13][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[13][0] + shift) / sigma)**2.)\
        + P[14][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - P[14][0] + shift) / sigma)**2.)))

    def S_fun(x, mS, shift, sigma):
        return (abs(mS) * \
        (+ S[0][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - S[0][0] + shift) / sigma)**2.)\
        + S[1][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - S[1][0] + shift) / sigma)**2.)\
        + S[2][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - S[2][0] + shift) / sigma)**2.)\
        + S[3][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - S[3][0] + shift) / sigma)**2.)\
        + S[4][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - S[4][0] + shift) / sigma)**2.)\
        + S[5][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - S[5][0] + shift) / sigma)**2.)\
        + S[6][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - S[6][0] + shift) / sigma)**2.)))
    
    def Z_fun(x, mZ, shift, sigma):
        return (abs(mZ) * \
        (+ Z[0][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[0][0] + shift) / sigma)**2.)\
        + Z[1][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[1][0] + shift) / sigma)**2.)\
        + Z[2][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[2][0] + shift) / sigma)**2.)\
        + Z[3][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[3][0] + shift) / sigma)**2.)\
        + Z[4][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[4][0] + shift) / sigma)**2.)\
        + Z[5][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[5][0] + shift) / sigma)**2.)\
        + Z[6][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[6][0] + shift) / sigma)**2.)\
        + Z[7][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[7][0] + shift) / sigma)**2.)\
        + Z[8][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[8][0] + shift) / sigma)**2.)\
        + Z[9][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[9][0] + shift) / sigma)**2.)\
        + Z[10][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[10][0] + shift) / sigma)**2.)\
        + Z[11][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[11][0] + shift) / sigma)**2.)\
        + Z[12][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[12][0] + shift) / sigma)**2.)\
        + Z[13][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[13][0] + shift) / sigma)**2.)\
        + Z[14][1]* 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - Z[14][0] + shift) / sigma)**2.)))

    def power_law2(x, amplitude, alpha):
        return amplitude * (x**-alpha)
    

    inter = interpolate.interp1d(wl_OII_cont, fl_OII_cont, kind='slinear')


    flux_F.append(integrate.trapz(F_fun(agn[:, 0], m.mF_0.value, m.shift_0.value, m.sigma_0.value)))
    flux_G.append(integrate.trapz(G_fun(agn[:, 0], m.mG_0.value, m.shift_0.value, m.sigma_0.value)))
    flux_P.append(integrate.trapz(P_fun(agn[:, 0], m.mP_0.value, m.shift_0.value, m.sigma_0.value)))
    flux_S.append(integrate.trapz(S_fun(agn[:, 0], m.mS_0.value, m.shift_0.value, m.sigma_0.value)))
    flux_Z.append(integrate.trapz(Z_fun(agn[:, 0], m.mZ_0.value, m.shift_0.value, m.sigma_0.value)))
    sigma_fe.append(m.sigma_0.value)
    shift_fe.append(m.shift_0.value)

    o = open('calculation/Fe_II_flux.csv', 'a')
    #for i in range(len(flux_F)):
    o.write(str(flux_F[q]) + '\t' + str(flux_G[q]) + '\t' + str(flux_P[q])\
     + '\t' + str(flux_S[q]) + '\t' + str(flux_Z[q]) + '\t' + str(sigma_fe[q])\
     + '\t' + str(shift_fe[q])+ '\n')
    
    o.close()



    plt.figure(j)
    plt.plot(agn[:, 0], agn[:, 1], 'b-', label='Observed Spectra')
    plt.plot(wl_hbeta, m(wl_hbeta), 'r-',label='Model')
    plt.plot(wl_halpha, m2(wl_halpha), 'r-')
    
    #--plt.plot(wl_OII, m3(wl_OII), 'r-')
    #--plt.plot(wl_OII, interpolate.spline(wl_OII_cont, fl_OII_cont, wl_OII, order=1), 'r-')
    
    plt.plot(wl_OII, inter(wl_OII), 'r-')

    plt.plot(wl_hbeta, power_law2(wl_hbeta, m.amplitude_1.value, m.alpha_1.value), 'g-', label='Power Law')
    plt.plot(wl_halpha, power_law2(wl_halpha, m2.amplitude_1.value, m2.alpha_1.value), 'g-')
    #--plt.plot(wl_OII, power_law2(wl_OII, m3.amplitude_1.value, m3.alpha_1.value), 'g-')



    plt.plot(wl_hbeta, F_fun(wl_hbeta, m.mF_0.value, m.shift_0.value, m.sigma_0.value), 'c-')
    plt.plot(wl_hbeta, G_fun(wl_hbeta, m.mG_0.value, m.shift_0.value, m.sigma_0.value), 'g-')
    plt.plot(wl_hbeta, P_fun(wl_hbeta, m.mP_0.value, m.shift_0.value, m.sigma_0.value), 'm-')
    plt.plot(wl_hbeta, S_fun(wl_hbeta, m.mS_0.value, m.shift_0.value, m.sigma_0.value), 'y-')
    plt.plot(wl_hbeta, Z_fun(wl_hbeta, m.mZ_0.value, m.shift_0.value, m.sigma_0.value), 'k-')
    
    plt.plot(wl_halpha, F_fun(wl_halpha, m2.mF_0.value, m2.shift_0.value, m2.sigma_0.value), 'c-')
    plt.plot(wl_halpha, G_fun(wl_halpha, m2.mG_0.value, m2.shift_0.value, m2.sigma_0.value), 'g-')
    plt.plot(wl_halpha, P_fun(wl_halpha, m2.mP_0.value, m2.shift_0.value, m2.sigma_0.value), 'm-')
    plt.plot(wl_halpha, S_fun(wl_halpha, m2.mS_0.value, m2.shift_0.value, m2.sigma_0.value), 'y-')
    plt.plot(wl_halpha, Z_fun(wl_halpha, m2.mZ_0.value, m2.shift_0.value, m2.sigma_0.value), 'k-')

    #--plt.plot(wl_OII, F_fun(wl_OII, m3.mF_0.value, m3.shift_0.value, m3.sigma_0.value), 'c-')
    #--plt.plot(wl_OII, G_fun(wl_OII, m3.mG_0.value, m3.shift_0.value, m3.sigma_0.value), 'g-')
    #--plt.plot(wl_OII, P_fun(wl_OII, m3.mP_0.value, m3.shift_0.value, m3.sigma_0.value), 'm-')
    #--plt.plot(wl_OII, S_fun(wl_OII, m3.mS_0.value, m3.shift_0.value, m3.sigma_0.value), 'y-')
    #--plt.plot(wl_OII, Z_fun(wl_OII, m3.mZ_0.value, m3.shift_0.value, m3.sigma_0.value), 'k-')


    plt.legend(loc='best', fontsize='small')
    
    plt.xlabel('Wavelength ($\\rm \AA$)', fontsize='x-large')
    plt.ylabel('Flux Density ($\\rm \\times 10^{-17} \ erg \ cm^{-2} s^{-1} \AA^{-1}$)', fontsize='x-large')
    plt.title('Continuum Subtraction\nspec-' + plate[j] + '-' + mjd[j] + '-' + fiberid[j] + '.fits')
    plt.savefig('figure2/fig2_'+ str(j) + '_spec-' + plate[j] + '-' + mjd[j] + '-' + fiberid[j])
    plt.close('all')
#'''========== Continuum Substraction =========='''

    fl_OII_sub = []
    fl_hbeta_sub = []
    fl_halpha_sub = []
    
    fl_OII_norm = []
    fl_hbeta_norm = []
    fl_halpha_norm = []


    for i in range(len(wl_OII)):
        fl_OII_sub.append(fl_OII[i] - inter(wl_OII)[i])
        fl_OII_norm.append(fl_OII[i] / inter(wl_OII)[i])
        
    for i in range(len(wl_hbeta)):
        fl_hbeta_sub.append(fl_hbeta[i] - m(wl_hbeta)[i])
        fl_hbeta_norm.append(fl_hbeta[i] / m(wl_hbeta)[i])

    for i in range(len(wl_halpha)):
        fl_halpha_sub.append(fl_halpha[i] - m(wl_halpha)[i])
        fl_halpha_norm.append(fl_halpha[i] / m(wl_halpha)[i])
        
    oii = 0
    hbeta = 0
    halpha = 0
    
    o = open('cont_sub/sub_'+'spec-' + plate[j] + '-' + mjd[j] + '-' + fiberid[j] + '.fits', 'w')
    o.write('#Wavelength (Angstroms)\Flux')
    o.write('\n\n\n')
    
    
    for i in range(len(agn)):
        
        if agn[i][0] >= min(wl_OII) and agn[i][0] <= max(wl_OII):
            o.write(str(agn[i][0]) + '\t' + str(fl_OII_sub[oii]) + '\n')
            oii += 1
            
        elif agn[i][0] >= min(wl_hbeta) and agn[i][0] <= max(wl_hbeta):
            o.write(str(agn[i][0]) + '\t' + str(fl_hbeta_sub[hbeta]) + '\n')
            hbeta += 1
        
        elif agn[i][0] >= min(wl_halpha) and agn[i][0] <= max(wl_halpha):
            o.write(str(agn[i][0]) + '\t' + str(fl_halpha_sub[halpha]) + '\n')
            halpha+=1
    
        else:
            o.write(str(agn[i][0]) + '\t' + str(0) + '\n')
            
    o.close()    
    #################
    oii = 0
    hbeta = 0
    halpha = 0
    
    o = open('cont_norm/norm_'+'spec-' + plate[j] + '-' + mjd[j] + '-' + fiberid[j] + '.fits', 'w')
    o.write('#Wavelength (Angstroms)\Flux')
    o.write('\n\n\n')
    
    
    for i in range(len(agn)):
        
        if agn[i][0] >= min(wl_OII) and agn[i][0] <= max(wl_OII):
            o.write(str(agn[i][0]) + '\t' + str(fl_OII_norm[oii]) + '\n')
            oii += 1
            
        elif agn[i][0] >= min(wl_hbeta) and agn[i][0] <= max(wl_hbeta):
            o.write(str(agn[i][0]) + '\t' + str(fl_hbeta_norm[hbeta]) + '\n')
            hbeta += 1
        
        elif agn[i][0] >= min(wl_halpha) and agn[i][0] <= max(wl_halpha):
            o.write(str(agn[i][0]) + '\t' + str(fl_halpha_norm[halpha]) + '\n')
            halpha+=1
    
        else:
            o.write(str(agn[i][0]) + '\t' + str(1) + '\n')# penting
            
    o.close()    
    
    fl_cont.append(np.mean(power_law2(wl_cont, m.amplitude_1.value, m.alpha_1.value)))
    
    f = open('calculation/continuum_flux_1.csv', 'a')
    #for i in range(len(fl_cont)):
    f.write(str(fl_cont[q]) + '\n')

    f.close()
    
    print'j = ', j
    
    q+=1
