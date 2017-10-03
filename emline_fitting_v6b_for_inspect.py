# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 11:32:12 2017

@author: Irham
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 09:43:58 2016

@author: irham
"""

import numpy as np
import pylab as plt
from astropy.modeling import fitting
from astropy.modeling.models import custom_model
from scipy import integrate
from scipy.stats import chisquare, f
#from PyAstronomy.pyasl import ftest



c = 299792458 * 10**-3 #km/s

o = open('calculation/m_class.csv', 'w')
o.write('m_class\n')
o.close()

o = open('calculation/m_class_alpha.csv', 'w')
o.write('m_class_alpha\n')
o.close()


o = open('calculation/line_flux.csv', 'w')
o.write('flux_[OII]_3727\tflux_n_H_beta_1\tflux_n_H_beta_2\
\tflux_b_H_beta_1\tflux_[OIII]_4959_1\tflux_[OIII]_4959_2\
\tflux_[OIII]_5007_1\tflux_[OIII]_5007_2\
\tflux_[OIII]_4959_wing\tflux_[OIII]_5007_wing\
\tflux_n_H_alpha_1\tflux_n_H_alpha_2\tflux_b_H_alpha_1\
\tflux_[NII]_6548_1\tflux_[NII]_6548_2\
\tflux_[NII]_6584_1\tflux_[NII]_6584_2\
\tflux_[SII]_6717_1\tflux_[SII]_6717_2\
\tflux_[SII]_6731_1\tflux_[SII]_6731_2\n')
o.close()


o = open('calculation/line_fwhm.csv', 'w')
o.write('sigma_[OII]_3727\tsigma_n_H_beta_1\tsigma_n_H_beta_2\
\tsigma_b_H_beta_1\tsigma_[OIII]_4959_1\tsigma_[OIII]_4959_2\
\tsigma_[OIII]_5007_1\tsigma_[OIII]_5007_2\
\tsigma_[OIII]_4959_wing\tsigma_[OIII]_5007_wing\
\tsigma_n_H_alpha_1\tsigma_n_H_alpha_2\tsigma_b_H_alpha_1\
\tsigma_[NII]_6548_1\tsigma_[NII]_6548_2\
\tsigma_[NII]_6584_1\tsigma_[NII]_6584_2\
\tsigma_[SII]_6717_1\tsigma_[SII]_6717_2\
\tsigma_[SII]_6731_1\tsigma_[SII]_6731_2\n')
o.close()

o = open('calculation/line_shift.csv', 'w')
o.write('shift_[OII]_3727\tshift_n_H_beta_1\tshift_n_H_beta_2\
\tshift_b_H_beta_1\tshift_[OIII]_4959_1\tshift_[OIII]_4959_2\
\tshift_[OIII]_5007_1\tshift_[OIII]_5007_2\
\tshift_[OIII]_4959_wing\tshift_[OIII]_5007_wing\
\tshift_n_H_alpha_1\tshift_n_H_alpha_2\tshift_b_H_alpha_1\
\tshift_[NII]_6548_1\tshift_[NII]_6548_2\
\tshift_[NII]_6584_1\tshift_[NII]_6584_2\
\tshift_[SII]_6717_1\tshift_[SII]_6717_2\
\tshift_[SII]_6731_1\tshift_[SII]_6731_2\n')
o.close()


o = open('calculation/line_amplitude.csv', 'w')
o.write('amp_n_H_alpha_1\tamp_n_H_alpha_2\tamp_b_H_alpha_1\n')
o.close()


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

#######################
sigma_OII_3727 = []

sigma_nHbeta_1 = []
sigma_nHbeta_2 = []
sigma_bHbeta_1 = []
sigma_OIII_4959_1 = []
sigma_OIII_5007_1 = []
sigma_OIII_4959_2 = []
sigma_OIII_5007_2 = []
sigma_OIII_4959_wing = []
sigma_OIII_5007_wing = []

sigma_nHalpha_1 = []
sigma_nHalpha_2 = []
sigma_bHalpha_1 = []
sigma_NII_6548_1 = []
sigma_NII_6584_1 = []
sigma_SII_6717_1 = []
sigma_SII_6731_1 = []
sigma_NII_6548_2 = []
sigma_NII_6584_2 = []
sigma_SII_6717_2 = []
sigma_SII_6731_2 = []

######################
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


def F_test (chisq1, chisq2, df1, df2):
    return ((chisq1-chisq2)/(df1-df2))/(chisq2/df2)

q = 0
M_type = 1

sample = np.loadtxt('MAIN_DATA.csv', skiprows=1, dtype=str, delimiter=',')

plate = []
mjd = []
fiberid = []
z = []
e_bv = []

skip = np.loadtxt('notes/skip.csv', delimiter=',', skiprows=1)


for j in range(35):#len(sample)):
    
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

    
    if j < 34:
        print 'skipped, j = ', j
        continue 
        
#    if j in skip[:, 1] or j == skip[0, 4]: #or j < 3455:#for until redshift 3.5
#        o = open('calculation/line_flux.csv', 'a')
#
#        o.write(str(0) + '\t' + str(0) + '\t' + str(0)\
#        + '\t' + str(0) + '\t' + str(0) + '\t' + str(0)\
#        + '\t' + str(0) + '\t' + str(0)\
#        + '\t' + str(0) + '\t' + str(0)\
#        + '\t' + str(0) + '\t' + str(0) + '\t' + str(0)\
#        + '\t' + str(0) + '\t' + str(0)\
#        + '\t' + str(0) + '\t' + str(0)\
#        + '\t' + str(0) + '\t' + str(0)\
#        + '\t' + str(0) + '\t' + str(0) + '\n')
#
#        o.close()
#
#        o = open('calculation/line_fwhm.csv', 'a')
#
#        o.write(str(0) + '\t' + str(0) + '\t' + str(0)\
#        + '\t' + str(0) + '\t' + str(0) + '\t' + str(0)\
#        + '\t' + str(0) + '\t' + str(0)\
#        + '\t' + str(0) + '\t' + str(0)\
#        + '\t' + str(0) + '\t' + str(0) + '\t' + str(0)\
#        + '\t' + str(0) + '\t' + str(0)\
#        + '\t' + str(0) + '\t' + str(0)\
#        + '\t' + str(0) + '\t' + str(0)\
#        + '\t' + str(0) + '\t' + str(0) + '\n')
#
#        o.close()
#
#        o = open('calculation/line_shift.csv', 'a')
#
#        o.write(str(0) + '\t' + str(0) + '\t' + str(0)\
#        + '\t' + str(0) + '\t' + str(0) + '\t' + str(0)\
#        + '\t' + str(0) + '\t' + str(0)\
#        + '\t' + str(0) + '\t' + str(0)\
#        + '\t' + str(0) + '\t' + str(0) + '\t' + str(0)\
#        + '\t' + str(0) + '\t' + str(0)\
#        + '\t' + str(0) + '\t' + str(0)\
#        + '\t' + str(0) + '\t' + str(0)\
#        + '\t' + str(0) + '\t' + str(0) + '\n')
#
#        o.close()
#
#        o = open('calculation/line_amplitude.csv', 'a')
#
#        o.write(str(0) + '\t' + str(0) + '\t' + str(0)\
#        + '\n')
#
#        o.close()
#
#        o = open('calculation/m_class.csv', 'a')
#        o.write(str(-999) + '\n')
#        o.close()
#
#        o = open('calculation/m_class_alpha.csv', 'a')
#        o.write(str(-999) + '\n')
#        o.close()
#
#
#        print 'skipped, j = ', j
#        continue


    agn = np.loadtxt('cont_sub/sub_spec-' + plate[j] + '-' + mjd[j] + '-' + fiberid[j] + '.fits')
    agn2 = np.loadtxt('cont_norm/norm_spec-' + plate[j] + '-' + mjd[j] + '-' + fiberid[j] + '.fits')
    
    wl_OII_3727 = []
    fl_OII_3727 = []

    wl_hbeta = []
    fl_hbeta = []

    wl_halpha = []
    fl_halpha = []

    wl_OIII_4959 = []
    fl_OIII_4959 = []

    wl_OIII_5007 = []
    fl_OIII_5007 = []
    
    
    wl_OII_3727_norm = []
    fl_OII_3727_norm = []

    wl_hbeta_norm = []
    fl_hbeta_norm = []

    wl_halpha_norm = []
    fl_halpha_norm = []

    wl_OIII_4959_norm = []
    fl_OIII_4959_norm = []

    wl_OIII_5007_norm = []
    fl_OIII_5007_norm = []    
    
    
    

    for i in range(len(agn)):#Hbeta window
        if agn[i][0] >= 3717 and agn[i][0] <=3740:
            wl_OII_3727.append(agn[i][0])
            fl_OII_3727.append(agn[i][1])
            
            wl_OII_3727_norm.append(agn2[i][0])
            fl_OII_3727_norm.append(agn2[i][1])
        
        elif agn[i][0] >= 6400 and agn[i][0] <=6800:
            wl_halpha.append(agn[i][0])
            fl_halpha.append(agn[i][1])
            
            wl_halpha_norm.append(agn2[i][0])
            fl_halpha_norm.append(agn2[i][1])
    
        elif agn[i][0] >= 4700 and agn[i][0] <=5100:#3747:
            wl_hbeta.append(agn[i][0])
            fl_hbeta.append(agn[i][1])   
            
            wl_hbeta_norm.append(agn2[i][0])
            fl_hbeta_norm.append(agn2[i][1])   
        else:
            continue
    
    fl_OII_3727_norm = np.array(fl_OII_3727_norm)
    fl_halpha_norm = np.array(fl_halpha_norm)
    fl_hbeta_norm = np.array(fl_hbeta_norm)

    
    @custom_model
    def broad_hbeta(x, amplitude_b1=np.mean(fl_hbeta), sigma_b1=3., shift_b=0):
        return ((abs(amplitude_b1)* 1./(abs(sigma_b1)*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - 4860 + shift_b) / abs(sigma_b1))**2.)))

    #narrow gaussian 1 hbeta
    @custom_model
    def narrow_line_1_hbeta(x, amplitude_n1=np.mean(fl_hbeta), amplitude_n2=np.mean(fl_hbeta), shift_n=0., sigma_n=1.):
        return (abs(amplitude_n1)* 1./(abs(sigma_n)*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - 4959 + shift_n) / abs(sigma_n))**2.)\
        + abs(3.*amplitude_n1)* 1./(abs(sigma_n)*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - 5007 + shift_n) / abs(sigma_n))**2.)\
        + abs(amplitude_n2)* 1./(abs(sigma_n)*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - 4860 + shift_n) / abs(sigma_n))**2.))

    #narrow gaussian 2 hbeta
    @custom_model
    def narrow_line_2_hbeta(x, amplitude_n1=np.mean(fl_hbeta), amplitude_n2=np.mean(fl_hbeta), shift_n=-5., sigma_n=1.):
        return (abs(amplitude_n1)* 1./(abs(sigma_n)*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - 4959 + shift_n) / abs(sigma_n))**2.)\
        + abs(3.*amplitude_n1)* 1./(abs(sigma_n)*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - 5007 + shift_n) / abs(sigma_n))**2.)\
        + abs(amplitude_n2)* 1./(abs(sigma_n)*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - 4860 + shift_n) / abs(sigma_n))**2.))


    #narrow gaussian 1 halpha
    @custom_model
    def narrow_line_1_halpha(x, amplitude_n3=np.mean(fl_halpha), amplitude_n4=np.mean(fl_halpha),\
    amplitude_n5=np.mean(fl_halpha), amplitude_n6=np.mean(fl_halpha),\
    shift_n=0., sigma_n=1.):
        return (abs(amplitude_n3)* 1./(abs(sigma_n)*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - 6563 + shift_n) / abs(sigma_n))**2.)\
        + abs(amplitude_n4)* 1./(abs(sigma_n)*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - 6548 + shift_n) / abs(sigma_n))**2.)\
        + abs((2.96)*amplitude_n4)* 1./(abs(sigma_n)*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - 6584 + shift_n) / abs(sigma_n))**2.)\

        + abs(amplitude_n5)* 1./(abs(sigma_n)*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - 6717 + shift_n) / abs(sigma_n))**2.)\
        + abs(amplitude_n6)* 1./(abs(sigma_n)*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - 6731 + shift_n) / abs(sigma_n))**2.))

    #narrow gaussian 2 halpha
    @custom_model
    def narrow_line_2_halpha(x, amplitude_n3=np.mean(fl_halpha), amplitude_n4=np.mean(fl_halpha),\
    amplitude_n5=np.mean(fl_halpha), amplitude_n6=np.mean(fl_halpha),\
    shift_n=-5., sigma_n=1.):
        return (abs(amplitude_n3)* 1./(abs(sigma_n)*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - 6563 + shift_n) / abs(sigma_n))**2.)\
        + abs(amplitude_n4)* 1./(abs(sigma_n)*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - 6548 + shift_n) / abs(sigma_n))**2.)\
        + abs((2.96)*amplitude_n4)* 1./(abs(sigma_n)*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - 6584 + shift_n) / abs(sigma_n))**2.)\

        + abs(amplitude_n5)* 1./(abs(sigma_n)*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - 6717 + shift_n) / abs(sigma_n))**2.)\
        + abs(amplitude_n6)* 1./(abs(sigma_n)*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - 6731 + shift_n) / abs(sigma_n))**2.))

    @custom_model
    def broad_halpha(x, amplitude_b1=np.mean(fl_hbeta), sigma_b1=3., shift_b=0):
        return ((abs(amplitude_b1)* 1./(abs(sigma_b1)*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - 6563 + shift_b) / abs(sigma_b1))**2.)))

    #[O III] wing
    @custom_model
    def oiii_wing(x, amplitude_n1=np.mean(fl_hbeta), shift_n=0., sigma_n=1.5):
        return (abs(amplitude_n1)* 1./(abs(sigma_n)*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - 4959 + shift_n) / abs(sigma_n))**2.)\
        + abs(3.*amplitude_n1)* 1./(abs(sigma_n)*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - 5007 + shift_n) / abs(sigma_n))**2.))


    @custom_model
    def n_oii(x, amplitude=np.mean(fl_OII_3727), sigma=1., shift=0.):
        return (abs(amplitude)* 1./(abs(sigma)*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - 3727 + shift) / abs(sigma))**2.))

    @custom_model
    def null(x):
        return 0

    Chi2_halpha = []
    dof_halpha = []
    F_pvalue_halpha = []    
    
    M_type = 0
    M_type == 1#narrow gaussian
    m_init3 = narrow_line_1_halpha() + narrow_line_2_halpha(amplitude_n3=0, amplitude_n4=0,\
    amplitude_n5=0, amplitude_n6=0)\
    + null() + broad_halpha(amplitude_b1=0.)
    m_init3.amplitude_n3_1.fixed = True #gaussian 2_c
    m_init3.amplitude_n4_1.fixed = True #gaussian 2_d
    m_init3.amplitude_n5_1.fixed = True #gaussian 2_e
    m_init3.amplitude_n6_1.fixed = True #gaussian 2_f
    m_init3.amplitude_b1_3.fixed = True #broad balmer
    fit = fitting.LevMarLSQFitter()
    m3_a = fit(m_init3, wl_halpha, fl_halpha, maxiter=10**3)        
    Chi2_halpha.append(chisquare(fl_halpha, m3_a(wl_halpha))[0])
    dof_halpha.append(len(wl_halpha) - 6)

    M_type == 2#narrow gaussian + narrow gaussian
    m_init3 = narrow_line_1_halpha() + narrow_line_2_halpha()\
    + null() + broad_halpha(amplitude_b1=0.)
    m_init3.amplitude_b1_3.fixed = True #broad balmer
    fit = fitting.LevMarLSQFitter()
    m3_b = fit(m_init3, wl_halpha, fl_halpha, maxiter=10**3)
    Chi2_halpha.append(chisquare(fl_halpha, m3_b(wl_halpha))[0])
    dof_halpha.append(len(wl_halpha) - 12)
    
    #print F_test(Chi2[0], Chi2[1], dof[0], dof[1])    
#    F_pvalue.append(ftest(Chi2[0], Chi2[1], dof[0], dof[1])['p-value']) #verified
#    print f.sf(F[0], dof[0], dof[1], loc=0, scale=1)
    
    M_type == 4#narrow gaussian + broad balmer
    m_init3 = narrow_line_1_halpha() + narrow_line_2_halpha(amplitude_n3=0, amplitude_n4=0,\
    amplitude_n5=0, amplitude_n6=0)\
    + null() + broad_halpha()
    m_init3.amplitude_n3_1.fixed = True #gaussian 2_c
    m_init3.amplitude_n4_1.fixed = True #gaussian 2_d
    m_init3.amplitude_n5_1.fixed = True #gaussian 2_e
    m_init3.amplitude_n6_1.fixed = True #gaussian 2_f
    fit = fitting.LevMarLSQFitter()
    m3_c = fit(m_init3, wl_halpha, fl_halpha, maxiter=10**3)        
    Chi2_halpha.append(chisquare(fl_halpha, m3_c(wl_halpha))[0])
    dof_halpha.append(len(wl_halpha) - 9)

    M_type == 6#narrow gaussian + narrow gaussian + broad balmer
    m_init3 = narrow_line_1_halpha() + narrow_line_2_halpha()\
    + null() + broad_halpha()
    fit = fitting.LevMarLSQFitter()
    m3_d = fit(m_init3, wl_halpha, fl_halpha, maxiter=10**3)        
    Chi2_halpha.append(chisquare(fl_halpha, m3_d(wl_halpha))[0])
    dof_halpha.append(len(wl_halpha) - 15)

#    F_pvalue_halpha.append(f.sf(F_test(Chi2_halpha[0], Chi2_halpha[2], dof_halpha[0], dof_halpha[2]), dof_halpha[0], dof_halpha[2]))
#    F_pvalue_halpha.append(f.sf(F_test(Chi2_halpha[2], Chi2_halpha[1], dof_halpha[2], dof_halpha[1]), dof_halpha[2], dof_halpha[1]))
#    F_pvalue_halpha.append(f.sf(F_test(Chi2_halpha[1], Chi2_halpha[3], dof_halpha[1], dof_halpha[3]), dof_halpha[1], dof_halpha[3]))

    F_pvalue_halpha.append(f.sf(F_test(Chi2_halpha[0], Chi2_halpha[1], dof_halpha[0], dof_halpha[1]), dof_halpha[0], dof_halpha[1]))
    F_pvalue_halpha.append(f.sf(F_test(Chi2_halpha[0], Chi2_halpha[2], dof_halpha[0], dof_halpha[2]), dof_halpha[0], dof_halpha[2]))
    F_pvalue_halpha.append(f.sf(F_test(Chi2_halpha[0], Chi2_halpha[3], dof_halpha[0], dof_halpha[3]), dof_halpha[0], dof_halpha[3]))
    F_pvalue_halpha.append(f.sf(F_test(Chi2_halpha[1], Chi2_halpha[0], dof_halpha[1], dof_halpha[0]), dof_halpha[1], dof_halpha[0]))
    F_pvalue_halpha.append(f.sf(F_test(Chi2_halpha[1], Chi2_halpha[2], dof_halpha[1], dof_halpha[2]), dof_halpha[1], dof_halpha[2]))    
    F_pvalue_halpha.append(f.sf(F_test(Chi2_halpha[1], Chi2_halpha[3], dof_halpha[1], dof_halpha[3]), dof_halpha[1], dof_halpha[3]))    
    F_pvalue_halpha.append(f.sf(F_test(Chi2_halpha[2], Chi2_halpha[0], dof_halpha[2], dof_halpha[0]), dof_halpha[2], dof_halpha[0]))        
    F_pvalue_halpha.append(f.sf(F_test(Chi2_halpha[2], Chi2_halpha[1], dof_halpha[2], dof_halpha[1]), dof_halpha[2], dof_halpha[1]))        
    F_pvalue_halpha.append(f.sf(F_test(Chi2_halpha[2], Chi2_halpha[3], dof_halpha[2], dof_halpha[3]), dof_halpha[2], dof_halpha[3]))            
    Chi2_red_halpha = np.array(Chi2_halpha)/np.array(dof_halpha)

    alpha_ = 1 - 99.73/100.

    def gauss(x, x0, amplitude, sigma, shift):
        return (abs(amplitude)* 1./(abs(sigma)*np.sqrt(2.*np.pi)) * np.exp(-0.5 * ((x - x0 + shift) / abs(sigma))**2.))

#    if F_pvalue_halpha[0] > alpha_\
#    and F_pvalue_halpha[1] > alpha_\
#    and F_pvalue_halpha[2] > alpha_:
#        m3 = m3_a
    
#    elif F_pvalue_halpha[3] > alpha_\
#    and F_pvalue_halpha[4] > alpha_\
#    and F_pvalue_halpha[5] > alpha_:
#        m3 = m3_b    

#    elif F_pvalue_halpha[6] > alpha_\
#    and F_pvalue_halpha[7] > alpha_\
#    and F_pvalue_halpha[8] > alpha_:
#        m3 = m3_c    
        
#    else:
#        m3 = m3_d
    wl_halpha = np.array(wl_halpha)

    if abs(m3_d.sigma_b1_3.value) * 2.3548 * c/(4860. - m3_d.shift_b_3.value) > 800.\
    and 6450. < 6563. - m3_d.shift_b_3.value < 6650.\
    and abs(m3_d.sigma_b1_3.value) > min(abs(m3_d.sigma_n_0.value), abs(m3_d.sigma_n_1.value))\
    \
    and 0.01 < integrate.trapz(gauss(wl_halpha, 6563., m3_d.amplitude_n3_0.value, m3_d.sigma_n_0.value, m3_d.shift_n_0.value))\
             / integrate.trapz(gauss(wl_halpha, 6563., m3_d.amplitude_n3_1.value, m3_d.sigma_n_1.value, m3_d.shift_n_1.value))\
        < 100.\
    \
    \
    and abs(c * ((m3_d.shift_n_0.value / (6563. - m3_d.shift_n_0.value))\
              -  (m3_d.shift_n_1.value / (6563. - m3_d.shift_n_1.value))))\
      / min(abs(m3_d.sigma_n_0.value)*2.3548*c/(6563. - m3_d.shift_n_0.value)\
      , abs(m3_d.sigma_n_1.value)* 2.3548*c/(6563. - m3_d.shift_n_1.value)) > 0.8\
    :
  
        m3 = m3_d

        o = open('calculation/m_class_alpha.csv', 'a')
        o.write(str(6) + '\n')
        o.close()


    elif abs(m3_c.sigma_b1_3.value) * 2.3548 * c/(4860. - m3_c.shift_b_3.value) > 800.\
    and 6450. < 6563. - m3_c.shift_b_3.value < 6650.\
    and abs(m3_c.sigma_b1_3.value) > min(abs(m3_c.sigma_n_0.value), abs(m3_c.sigma_n_1.value))\
    \
    :
  
        m3 = m3_c

        o = open('calculation/m_class_alpha.csv', 'a')
        o.write(str(4) + '\n')
        o.close()

    elif 0.01 < integrate.trapz(gauss(wl_halpha, 6563., m3_b.amplitude_n3_0.value, m3_b.sigma_n_0.value, m3_b.shift_n_0.value))\
             / integrate.trapz(gauss(wl_halpha, 6563., m3_b.amplitude_n3_1.value, m3_b.sigma_n_1.value, m3_b.shift_n_1.value))\
        < 100.\
    \
    \
    and abs(c * ((m3_b.shift_n_0.value / (6563. - m3_b.shift_n_0.value))\
              -  (m3_b.shift_n_1.value / (6563. - m3_b.shift_n_1.value))))\
      / min(abs(m3_b.sigma_n_0.value)*2.3548*c/(6563. - m3_b.shift_n_0.value)\
      , abs(m3_b.sigma_n_1.value)* 2.3548*c/(6563. - m3_b.shift_n_1.value)) > 0.8\
    :
  
        m3 = m3_b

        o = open('calculation/m_class_alpha.csv', 'a')
        o.write(str(2) + '\n')
        o.close()

    else:
        m3 = m3_a

        o = open('calculation/m_class_alpha.csv', 'a')
        o.write(str(1) + '\n')
        o.close()

#    m3 = m3_a
#    if F_pvalue_halpha[0] < alpha_:
#        m3 = m3_c
#    if F_pvalue_halpha[1] < alpha_:
#        m3 = m3_b
#    if F_pvalue_halpha[2] < alpha_:
#        m3 = m3_d


############ FOR HBETA ########################
    Chi2_hbeta = []
    dof_hbeta = []
    F_pvalue_hbeta = []    



    M_type = 8
    M_type == 1#narrow gaussian
    m_init = narrow_line_1_hbeta() + narrow_line_2_hbeta(amplitude_n1=0, amplitude_n2=0)\
    + oiii_wing(amplitude_n1=0.) + broad_hbeta(amplitude_b1=0.)
    m_init.amplitude_n1_1.fixed = True #gaussian 2_a
    m_init.amplitude_n2_1.fixed = True #gaussian 2_b
    m_init.amplitude_n1_2.fixed = True #oiii wing
    m_init.amplitude_b1_3.fixed = True #broad balmer
    fit = fitting.LevMarLSQFitter()
    m_a = fit(m_init, wl_hbeta, fl_hbeta, maxiter=10**3)        
    Chi2_hbeta.append(chisquare(fl_hbeta, m_a(wl_hbeta))[0])
    dof_hbeta.append(len(wl_hbeta) - 4)

    
    M_type == 2#narrow gaussian + narrow gaussian
    m_init = narrow_line_1_hbeta() + narrow_line_2_hbeta()\
    + oiii_wing(amplitude_n1=0.) + broad_hbeta(amplitude_b1=0.)
    m_init.amplitude_n1_2.fixed = True #oiii wing
    m_init.amplitude_b1_3.fixed = True #broad balmer
    fit = fitting.LevMarLSQFitter()
    m_b = fit(m_init, wl_hbeta, fl_hbeta, maxiter=10**3)
    Chi2_hbeta.append(chisquare(fl_hbeta, m_b(wl_hbeta))[0])
    dof_hbeta.append(len(wl_hbeta) - 8)


    M_type == 3#narrow gaussian + [O III] wings
    m_init = narrow_line_1_hbeta() + narrow_line_2_hbeta(amplitude_n1=0, amplitude_n2=0)\
    + oiii_wing() + broad_hbeta(amplitude_b1=0.)
    m_init.amplitude_n1_1.fixed = True #gaussian 2_a
    m_init.amplitude_n2_1.fixed = True #gaussian 2_b
    m_init.amplitude_b1_3.fixed = True #broad balmer
    fit = fitting.LevMarLSQFitter()
    m_c = fit(m_init, wl_hbeta, fl_hbeta, maxiter=10**3)
    Chi2_hbeta.append(chisquare(fl_hbeta, m_c(wl_hbeta))[0])
    dof_hbeta.append(len(wl_hbeta) - 7)

    M_type == 4#narrow gaussian + broad Balmer
    m_init = narrow_line_1_hbeta() + narrow_line_2_hbeta(amplitude_n1=0, amplitude_n2=0)\
    + oiii_wing(amplitude_n1=0.) + broad_hbeta()
    m_init.amplitude_n1_1.fixed = True #gaussian 2_a
    m_init.amplitude_n2_1.fixed = True #gaussian 2_b
    m_init.amplitude_n1_2.fixed = True #oiii wing
    fit = fitting.LevMarLSQFitter()
    m_d = fit(m_init, wl_hbeta, fl_hbeta, maxiter=10**3)        
    Chi2_hbeta.append(chisquare(fl_hbeta, m_d(wl_hbeta))[0])
    dof_hbeta.append(len(wl_hbeta) - 7)


    M_type == 5#narrow gaussian + [O III] wings + broad Balmer
    m_init = narrow_line_1_hbeta() + narrow_line_2_hbeta(amplitude_n1=0, amplitude_n2=0)\
    + oiii_wing() + broad_hbeta()
    m_init.amplitude_n1_1.fixed = True #gaussian 2_a
    m_init.amplitude_n2_1.fixed = True #gaussian 2_b
    fit = fitting.LevMarLSQFitter()
    m_e = fit(m_init, wl_hbeta, fl_hbeta, maxiter=10**3)        
    Chi2_hbeta.append(chisquare(fl_hbeta, m_e(wl_hbeta))[0])
    dof_hbeta.append(len(wl_hbeta) - 10)

    M_type == 6#narrow gaussian + narrow gaussian + broad Balmer
    m_init = narrow_line_1_hbeta() + narrow_line_2_hbeta()\
    + oiii_wing(amplitude_n1=0.) + broad_hbeta()
    m_init.amplitude_n1_2.fixed = True #oiii wing
    fit = fitting.LevMarLSQFitter()
    m_f = fit(m_init, wl_hbeta, fl_hbeta, maxiter=10**3)        
    Chi2_hbeta.append(chisquare(fl_hbeta, m_f(wl_hbeta))[0])
    dof_hbeta.append(len(wl_hbeta) - 11)


    M_type == 7#narrow gaussian + narrow gaussian + [O III] wings
    m_init = narrow_line_1_hbeta() + narrow_line_2_hbeta()\
    + oiii_wing() + broad_hbeta(amplitude_b1=0.)
    m_init.amplitude_b1_3.fixed = True #broad balmer        
    fit = fitting.LevMarLSQFitter()
    m_g = fit(m_init, wl_hbeta, fl_hbeta, maxiter=10**3)        
    Chi2_hbeta.append(chisquare(fl_hbeta, m_g(wl_hbeta))[0])
    dof_hbeta.append(len(wl_hbeta) - 11)


    M_type == 8#narrow gaussian + narrow gaussian + [O III] wings + broad Balmer
    m_init = narrow_line_1_hbeta() + narrow_line_2_hbeta()\
    + oiii_wing() + broad_hbeta()
    fit = fitting.LevMarLSQFitter()
    m_h = fit(m_init, wl_hbeta, fl_hbeta, maxiter=10**3)
    Chi2_hbeta.append(chisquare(fl_hbeta, m_h(wl_hbeta))[0])
    dof_hbeta.append(len(wl_hbeta) - 14)

    for k in range(8):
        for l in range(8):
            if k != l:
                F_pvalue_hbeta.append(f.sf(F_test(Chi2_hbeta[k], Chi2_hbeta[l], dof_hbeta[k], dof_hbeta[l]), dof_hbeta[k], dof_hbeta[l]))
                #print k+1, l+1, f.sf(F_test(Chi2_hbeta[k], Chi2_hbeta[l], dof_hbeta[k], dof_hbeta[l]), dof_hbeta[k], dof_hbeta[l])




    wl_hbeta = np.array(wl_hbeta)    

#    and abs(c * ((m_h.shift_n_0.value / (5007. - m_h.shift_n_0.value))\
#              -  (m_h.shift_n_1.value / (5007. - m_h.shift_n_1.value)))) > 200.\

#    and abs(3.*m_h.amplitude_n1_2.value) < abs(3.*m_h.amplitude_n1_0.value)\
#    and abs(3.*m_h.amplitude_n1_2.value) < abs(3.*m_h.amplitude_n1_1.value):
    
    #Let's begin
    if abs(m_h.sigma_b1_3.value) * 2.3548 * c/(4860. - m_h.shift_b_3.value) > 800.\
    and 4770. < 4860. - m_h.shift_b_3.value < 4950.\
    and abs(m_h.sigma_b1_3.value) > min(abs(m_h.sigma_n_0.value), abs(m_h.sigma_n_1.value))\
    \
    and 0.01 < integrate.trapz(gauss(wl_hbeta, 5007., 3.*m_h.amplitude_n1_0.value, m_h.sigma_n_0.value, m_h.shift_n_0.value))\
             / integrate.trapz(gauss(wl_hbeta, 5007., 3.*m_h.amplitude_n1_1.value, m_h.sigma_n_1.value, m_h.shift_n_1.value))\
        < 100.\
    \
    \
    and abs(c * ((m_h.shift_n_0.value / (5007. - m_h.shift_n_0.value))\
              -  (m_h.shift_n_1.value / (5007. - m_h.shift_n_1.value))))\
      / min(abs(m_h.sigma_n_0.value)*2.3548*c/(5007. - m_h.shift_n_0.value)\
      , abs(m_h.sigma_n_1.value)* 2.3548*c/(5007. - m_h.shift_n_1.value)) > 0.8\
    \
    and abs(m_h.shift_n_2.value < 20.)\
    and abs(m_h.sigma_n_2.value) > max(abs(m_h.sigma_n_0.value), abs(m_h.sigma_n_1.value))\
    :

  
        m = m_h

        o = open('calculation/m_class.csv', 'a')
        o.write(str(8) + '\n')
        o.close() 

    elif 0.01 < integrate.trapz(gauss(wl_hbeta, 5007., 3.*m_g.amplitude_n1_0.value, m_g.sigma_n_0.value, m_g.shift_n_0.value))\
              / integrate.trapz(gauss(wl_hbeta, 5007., 3.*m_g.amplitude_n1_1.value, m_g.sigma_n_1.value, m_g.shift_n_1.value))\
              < 100.\
    \
    \
    and abs(c * ((m_g.shift_n_0.value / (5007. - m_g.shift_n_0.value))\
              -  (m_g.shift_n_1.value / (5007. - m_g.shift_n_1.value))))\
      / min(abs(m_g.sigma_n_0.value)*2.3548*c/(5007. - m_g.shift_n_0.value)\
      , abs(m_g.sigma_n_1.value)* 2.3548*c/(5007. - m_g.shift_n_1.value)) > 0.8\
    \
    and abs(m_g.shift_n_2.value < 20.)\
    and abs(m_g.sigma_n_2.value) > max(abs(m_g.sigma_n_0.value), abs(m_g.sigma_n_1.value))\
    :
        m = m_g

        o = open('calculation/m_class.csv', 'a')
        o.write(str(7) + '\n')
        o.close() 

    elif abs(m_f.sigma_b1_3.value) * 2.3548 * c/(4860. - m_f.shift_b_3.value) > 800.\
    and 4770. < 4860. - m_f.shift_b_3.value < 4950.\
    and abs(m_f.sigma_b1_3.value) > min(abs(m_f.sigma_n_0.value), abs(m_f.sigma_n_1.value))\
    \
    and 0.01 < integrate.trapz(gauss(wl_hbeta, 5007., 3.*m_f.amplitude_n1_0.value, m_f.sigma_n_0.value, m_f.shift_n_0.value))\
             / integrate.trapz(gauss(wl_hbeta, 5007., 3.*m_f.amplitude_n1_1.value, m_f.sigma_n_1.value, m_f.shift_n_1.value))\
        < 100.\
    \
    \
    and abs(c * ((m_f.shift_n_0.value / (5007. - m_f.shift_n_0.value))\
              -  (m_f.shift_n_1.value / (5007. - m_f.shift_n_1.value))))\
      / min(abs(m_f.sigma_n_0.value)*2.3548*c/(5007. - m_f.shift_n_0.value)\
      , abs(m_f.sigma_n_1.value)* 2.3548*c/(5007. - m_f.shift_n_1.value)) > 0.8\
    :
  
        m = m_f

        o = open('calculation/m_class.csv', 'a')
        o.write(str(6) + '\n')
        o.close() 

    #Let's begin
    elif abs(m_e.sigma_b1_3.value) * 2.3548 * c/(4860. - m_e.shift_b_3.value) > 800.\
    and 4770. < 4860. - m_e.shift_b_3.value < 4950.\
    and abs(m_e.sigma_b1_3.value) > min(abs(m_e.sigma_n_0.value), abs(m_e.sigma_n_1.value))\
    \
    and abs(m_e.shift_n_2.value < 20.)\
    and abs(m_e.sigma_n_2.value) > max(abs(m_e.sigma_n_0.value), abs(m_e.sigma_n_1.value))\
    :
        m = m_e

        o = open('calculation/m_class.csv', 'a')
        o.write(str(5) + '\n')
        o.close()

    elif abs(m_d.sigma_b1_3.value) * 2.3548 * c/(4860. - m_d.shift_b_3.value) > 800.\
    and 4770. < 4860. - m_d.shift_b_3.value < 4950.\
    and abs(m_d.sigma_b1_3.value) > min(abs(m_d.sigma_n_0.value), abs(m_d.sigma_n_1.value))\
    :
        m = m_d

        o = open('calculation/m_class.csv', 'a')
        o.write(str(4) + '\n')
        o.close()

    elif abs(m_c.shift_n_2.value < 20.)\
    and abs(m_c.sigma_n_2.value) > max(abs(m_c.sigma_n_0.value), abs(m_c.sigma_n_1.value))\
    :
        m = m_c

        o = open('calculation/m_class.csv', 'a')
        o.write(str(3) + '\n')
        o.close()

    elif 0.01 < integrate.trapz(gauss(wl_hbeta, 5007., 3.*m_b.amplitude_n1_0.value, m_b.sigma_n_0.value, m_b.shift_n_0.value))\
             / integrate.trapz(gauss(wl_hbeta, 5007., 3.*m_b.amplitude_n1_1.value, m_b.sigma_n_1.value, m_b.shift_n_1.value))\
        < 100.\
    \
    \
    and abs(c * ((m_b.shift_n_0.value / (5007. - m_b.shift_n_0.value))\
              -  (m_b.shift_n_1.value / (5007. - m_b.shift_n_1.value))))\
      / min(abs(m_b.sigma_n_0.value)*2.3548*c/(5007. - m_b.shift_n_0.value)\
      , abs(m_b.sigma_n_1.value)* 2.3548*c/(5007. - m_b.shift_n_1.value)) > 0.8\
    :

        m = m_b

        o = open('calculation/m_class.csv', 'a')
        o.write(str(2) + '\n')
        o.close()


    else:
        m = m_a

        o = open('calculation/m_class.csv', 'a')
        o.write(str(1) + '\n')
        o.close()
    
    if j == 34:
        m = m_b
    if j == 35:
        m = m_a
    #m2 ke m7 harus ada F-test
    m_init2 = n_oii()
    m2 = fit(m_init2, wl_OII_3727, fl_OII_3727, maxiter=10**3)
    perr_oii = np.sqrt(np.diag(fit.fit_info['param_cov']))
    

        
    b_hbeta_1 = []
    
    n_hbeta_1 = []
    n_hbeta_2 = []
    n_oiii_4959_1 = []
    n_oiii_5007_1 = []
    n_oiii_4959_2 = []
    n_oiii_5007_2 = []

    oiii_4959_wing = []
    oiii_5007_wing = []

    n_oii_3727 = []
    
    b_halpha_1 = []
    
    n_halpha_1 = []
    n_nii_6548_1 = []
    n_nii_6584_1 = []
    n_sii_6717_1 = []
    n_sii_6731_1 = []

    n_halpha_2 = []
    n_nii_6548_2 = []
    n_nii_6584_2 = []
    n_sii_6717_2 = []
    n_sii_6731_2 = []    


    for i in range(len(wl_hbeta)):
        n_oiii_4959_1.append(gauss(wl_hbeta[i], 4959., m.amplitude_n1_0.value, m.sigma_n_0.value, m.shift_n_0.value))        
        n_oiii_5007_1.append(gauss(wl_hbeta[i], 5007., 3.*m.amplitude_n1_0.value, m.sigma_n_0.value, m.shift_n_0.value))        
        n_hbeta_1.append(gauss(wl_hbeta[i], 4860., m.amplitude_n2_0.value, m.sigma_n_0.value, m.shift_n_0.value))

        n_oiii_4959_2.append(gauss(wl_hbeta[i], 4959., m.amplitude_n1_1.value, m.sigma_n_1.value, m.shift_n_1.value))        
        n_oiii_5007_2.append(gauss(wl_hbeta[i], 5007., 3.*m.amplitude_n1_1.value, m.sigma_n_1.value, m.shift_n_1.value))        
        n_hbeta_2.append(gauss(wl_hbeta[i], 4860., m.amplitude_n2_1.value, m.sigma_n_1.value, m.shift_n_1.value))

        oiii_4959_wing.append(gauss(wl_hbeta[i], 4959., m.amplitude_n1_2.value, m.sigma_n_2.value, m.shift_n_2.value))        
        oiii_5007_wing.append(gauss(wl_hbeta[i], 5007., 3.*m.amplitude_n1_2.value, m.sigma_n_2.value, m.shift_n_2.value))        

    
        b_hbeta_1.append(gauss(wl_hbeta[i], 4860., m.amplitude_b1_3.value, m.sigma_b1_3.value, m.shift_b_3.value))    

    for i in range(len(wl_halpha)):
        #######
        b_halpha_1.append(gauss(wl_halpha[i], 6563., m3.amplitude_b1_3.value, m3.sigma_b1_3.value, m3.shift_b_3.value))
#        + gauss(wl_halpha[i], 6563., m3.amplitude_b2_1.value, m3.sigma_b2_1.value, m3.shift_b2_1.value))
        
        n_halpha_1.append(gauss(wl_halpha[i], 6563., m3.amplitude_n3_0.value, m3.sigma_n_0.value, m3.shift_n_0.value))

        n_nii_6548_1.append(gauss(wl_halpha[i], 6548, m3.amplitude_n4_0.value, m3.sigma_n_0.value, m3.shift_n_0.value))
        n_nii_6584_1.append(gauss(wl_halpha[i], 6584., (2.96)*m3.amplitude_n4_0.value, m3.sigma_n_0.value, m3.shift_n_0.value))
    
        n_sii_6717_1.append(gauss(wl_halpha[i], 6717., m3.amplitude_n5_0.value, m3.sigma_n_0.value, m3.shift_n_0.value))
        n_sii_6731_1.append(gauss(wl_halpha[i], 6731., m3.amplitude_n6_0.value, m3.sigma_n_0.value, m3.shift_n_0.value))

        n_halpha_2.append(gauss(wl_halpha[i], 6563., m3.amplitude_n3_1.value, m3.sigma_n_1.value, m3.shift_n_1.value))

        n_nii_6548_2.append(gauss(wl_halpha[i], 6548, m3.amplitude_n4_1.value, m3.sigma_n_1.value, m3.shift_n_1.value))
        n_nii_6584_2.append(gauss(wl_halpha[i], 6584., (2.96)*m3.amplitude_n4_1.value, m3.sigma_n_1.value, m3.shift_n_1.value))
    
        n_sii_6717_2.append(gauss(wl_halpha[i], 6717., m3.amplitude_n5_1.value, m3.sigma_n_1.value, m3.shift_n_1.value))
        n_sii_6731_2.append(gauss(wl_halpha[i], 6731., m3.amplitude_n6_1.value, m3.sigma_n_1.value, m3.shift_n_1.value))


    flux_OII_3727.append(integrate.trapz(m2(wl_OII_3727)))
 
    flux_nHbeta_1.append(integrate.trapz(n_hbeta_1))
    flux_nHbeta_2.append(integrate.trapz(n_hbeta_2))
    flux_bHbeta_1.append(integrate.trapz(b_hbeta_1))
    flux_OIII_4959_1.append(integrate.trapz(n_oiii_4959_1))
    flux_OIII_5007_1.append(integrate.trapz(n_oiii_5007_1))    
    flux_OIII_4959_2.append(integrate.trapz(n_oiii_4959_2))
    flux_OIII_5007_2.append(integrate.trapz(n_oiii_5007_2))    
    flux_OIII_4959_wing.append(integrate.trapz(oiii_4959_wing))
    flux_OIII_5007_wing.append(integrate.trapz(oiii_5007_wing))    

    flux_nHalpha_1.append(integrate.trapz(n_halpha_1))    
    flux_nHalpha_2.append(integrate.trapz(n_halpha_2))    
    flux_bHalpha_1.append(integrate.trapz(b_halpha_1))
    flux_NII_6548_1.append(integrate.trapz(n_nii_6548_1))
    flux_NII_6584_1.append(integrate.trapz(n_nii_6584_1))    
    flux_SII_6717_1.append(integrate.trapz(n_sii_6717_1))
    flux_SII_6731_1.append(integrate.trapz(n_sii_6731_1))
    flux_NII_6548_2.append(integrate.trapz(n_nii_6548_2))
    flux_NII_6584_2.append(integrate.trapz(n_nii_6584_2))    
    flux_SII_6717_2.append(integrate.trapz(n_sii_6717_2))
    flux_SII_6731_2.append(integrate.trapz(n_sii_6731_2))

#############################################################
    sigma_OII_3727.append(abs(m2.sigma.value))# * 2.3548 * c/(3727.+ m2.shift.value))
    sigma_OIII_4959_1.append(abs(m.sigma_n_0.value))
    sigma_OIII_5007_1.append(abs(m.sigma_n_0.value))        
    sigma_nHbeta_1.append(abs(m.sigma_n_0.value))
    sigma_OIII_4959_2.append(abs(m.sigma_n_1.value))        
    sigma_OIII_5007_2.append(abs(m.sigma_n_1.value))        
    sigma_nHbeta_2.append(abs(m.sigma_n_1.value))
    sigma_OIII_4959_wing.append(abs(m.sigma_n_2.value))        
    sigma_OIII_5007_wing.append(abs(m.sigma_n_2.value))            
    sigma_bHbeta_1.append(abs(m.sigma_b1_3.value))    

    sigma_bHalpha_1.append(abs(m3.sigma_b1_3.value))        
    sigma_nHalpha_1.append(abs(m3.sigma_n_0.value))
    sigma_NII_6548_1.append(abs(m3.sigma_n_0.value))
    sigma_NII_6584_1.append(abs(m3.sigma_n_0.value))    
    sigma_SII_6717_1.append(abs(m3.sigma_n_0.value))
    sigma_SII_6731_1.append(abs(m3.sigma_n_0.value))
    sigma_nHalpha_2.append(abs(m3.sigma_n_1.value))
    sigma_NII_6548_2.append(abs(m3.sigma_n_1.value))
    sigma_NII_6584_2.append(abs(m3.sigma_n_1.value))    
    sigma_SII_6717_2.append(abs(m3.sigma_n_1.value))
    sigma_SII_6731_2.append(abs(m3.sigma_n_1.value))

###########################################################
    shift_OII_3727.append(m2.shift.value)# * 2.3548 * c/(3727.+ m2.shift.value))
    shift_OIII_4959_1.append(m.shift_n_0.value)        
    shift_OIII_5007_1.append(m.shift_n_0.value)        
    shift_nHbeta_1.append(m.shift_n_0.value)
    shift_OIII_4959_2.append(m.shift_n_1.value)        
    shift_OIII_5007_2.append(m.shift_n_1.value)        
    shift_nHbeta_2.append(m.shift_n_1.value)
    shift_OIII_4959_wing.append(m.shift_n_2.value)        
    shift_OIII_5007_wing.append(m.shift_n_2.value)
    shift_bHbeta_1.append(m.shift_b_3.value)    

    shift_bHalpha_1.append(m3.shift_b_3.value)        
    shift_nHalpha_1.append(m3.shift_n_0.value)
    shift_NII_6548_1.append(m3.shift_n_0.value)
    shift_NII_6584_1.append(m3.shift_n_0.value)    
    shift_SII_6717_1.append(m3.shift_n_0.value)
    shift_SII_6731_1.append(m3.shift_n_0.value)
    shift_nHalpha_2.append(m3.shift_n_1.value)
    shift_NII_6548_2.append(m3.shift_n_1.value)
    shift_NII_6584_2.append(m3.shift_n_1.value)    
    shift_SII_6717_2.append(m3.shift_n_1.value)
    shift_SII_6731_2.append(m3.shift_n_1.value)

####################################################
    amp_bHalpha_1.append(max(b_halpha_1))
    amp_nHalpha_1.append(max(n_halpha_1))
    amp_nHalpha_2.append(max(n_halpha_2))

    o = open('calculation/line_flux.csv', 'a')


    #for i in range(len(flux_nHbeta)):
    o.write(str(flux_OII_3727[q]) + '\t' + str(flux_nHbeta_1[q]) + '\t' + str(flux_nHbeta_2[q])\
    + '\t' + str(flux_bHbeta_1[q]) + '\t' + str(flux_OIII_4959_1[q]) + '\t' + str(flux_OIII_4959_2[q])\
    + '\t' + str(flux_OIII_5007_1[q])+ '\t' + str(flux_OIII_5007_2[q])\
    + '\t' + str(flux_OIII_4959_wing[q])+ '\t' + str(flux_OIII_5007_wing[q])\
    + '\t' + str(flux_nHalpha_1[q])+ '\t' + str(flux_nHalpha_2[q]) + '\t' + str(flux_bHalpha_1[q])\
    + '\t' + str(flux_NII_6548_1[q]) + '\t' + str(flux_NII_6548_2[q])\
    + '\t' + str(flux_NII_6584_1[q])  + '\t' + str(flux_NII_6584_2[q])\
    + '\t' + str(flux_SII_6717_1[q]) + '\t' + str(flux_SII_6717_2[q])\
    + '\t' + str(flux_SII_6731_1[q]) + '\t' + str(flux_SII_6731_2[q]) + '\n')

    o.close()

    o = open('calculation/line_fwhm.csv', 'a')


    #for i in range(len(flux_nHbeta)):
    o.write(str(sigma_OII_3727[q]) + '\t' + str(sigma_nHbeta_1[q]) + '\t' + str(sigma_nHbeta_2[q])\
    + '\t' + str(sigma_bHbeta_1[q]) + '\t' + str(sigma_OIII_4959_1[q]) + '\t' + str(sigma_OIII_4959_2[q])\
    + '\t' + str(sigma_OIII_5007_1[q])+ '\t' + str(sigma_OIII_5007_2[q])\
    + '\t' + str(sigma_OIII_4959_wing[q])+ '\t' + str(sigma_OIII_5007_wing[q])\
    + '\t' + str(sigma_nHalpha_1[q])+ '\t' + str(sigma_nHalpha_2[q]) + '\t' + str(sigma_bHalpha_1[q])\
    + '\t' + str(sigma_NII_6548_1[q]) + '\t' + str(sigma_NII_6548_2[q])\
    + '\t' + str(sigma_NII_6584_1[q])  + '\t' + str(sigma_NII_6584_2[q])\
    + '\t' + str(sigma_SII_6717_1[q]) + '\t' + str(sigma_SII_6717_2[q])\
    + '\t' + str(sigma_SII_6731_1[q]) + '\t' + str(sigma_SII_6731_2[q]) + '\n')

    o.close()

    o = open('calculation/line_shift.csv', 'a')


    #for i in range(len(flux_nHbeta)):
    o.write(str(shift_OII_3727[q]) + '\t' + str(shift_nHbeta_1[q]) + '\t' + str(shift_nHbeta_2[q])\
    + '\t' + str(shift_bHbeta_1[q]) + '\t' + str(shift_OIII_4959_1[q]) + '\t' + str(shift_OIII_4959_2[q])\
    + '\t' + str(shift_OIII_5007_1[q])+ '\t' + str(shift_OIII_5007_2[q])\
    + '\t' + str(shift_OIII_4959_wing[q])+ '\t' + str(shift_OIII_5007_wing[q])\
    + '\t' + str(shift_nHalpha_1[q])+ '\t' + str(shift_nHalpha_2[q]) + '\t' + str(shift_bHalpha_1[q])\
    + '\t' + str(shift_NII_6548_1[q]) + '\t' + str(shift_NII_6548_2[q])\
    + '\t' + str(shift_NII_6584_1[q])  + '\t' + str(shift_NII_6584_2[q])\
    + '\t' + str(shift_SII_6717_1[q]) + '\t' + str(shift_SII_6717_2[q])\
    + '\t' + str(shift_SII_6731_1[q]) + '\t' + str(shift_SII_6731_2[q]) + '\n')

    o.close()


    o = open('calculation/line_amplitude.csv', 'a')

    #for i in range(len(flux_nHbeta)):
    o.write(str(amp_nHalpha_1[q])+ '\t' + str(amp_nHalpha_2[q]) + '\t' + str(amp_bHalpha_1[q]) + '\n')

    o.close()


    
#    plt.figure(j)
#    plt.plot(wl_hbeta, fl_hbeta, 'b-', linewidth=1.5, label='SDSS Spectrum')
##
##
##    plt.plot(wl_hbeta, n_oiii_4959_1, 'c-')
##    plt.plot(wl_hbeta, n_oiii_5007_1, 'c-')
##    plt.plot(wl_hbeta, n_hbeta_1, 'c-')
##
##    plt.plot(wl_hbeta, n_oiii_4959_2, 'g-')
##    plt.plot(wl_hbeta, n_oiii_5007_2, 'g-')
##    plt.plot(wl_hbeta, n_hbeta_2, 'g-')
##
##    plt.plot(wl_hbeta, oiii_4959_wing, 'y-')
##    plt.plot(wl_hbeta, oiii_5007_wing, 'y-')
##
##    plt.plot(wl_hbeta, b_hbeta_1, 'k-')
#    ap = []
#    for i in range(len(wl_hbeta)):
#        ap.append(m(wl_hbeta[i])+np.random.normal(loc=np.median(fl_hbeta), scale=1.5))
##    plt.plot(wl_hbeta, m(wl_hbeta), 'r-', label='MRES Spectrum')
#    plt.plot(wl_hbeta, ap, 'r-', label='MRES Spectrum')
#    plt.xlim(4800, 5050)
#
#    plt.xlabel('Wavelength ($\\rm \AA$)', fontsize='x-large')
#    plt.ylabel('Flux Density ($\\rm \\times 10^{-17} \ erg \ cm^{-2} s^{-1} \AA^{-1}$)', fontsize='x-large')
#    plt.title('Emission Lines Fitting\nspec-' + plate[j] + '-' + mjd[j] + '-' + fiberid[j] + '.fits')
#    plt.savefig('figure3/fig3_'+ str(j) +'_spec-' + plate[j] + '-' + mjd[j] + '-' + fiberid[j])
#    plt.legend(loc='best')
##    plt.close('all')
#    
#    plt.plot(wl_OII_3727, fl_OII_3727, 'b-', label='Observed Spectra')
#    plt.plot(wl_OII_3727, m2(wl_OII_3727), 'r-', label='Model')

##################################################

    plt.figure(j)
    plt.plot(wl_halpha, fl_halpha, 'b-', linewidth=1.5, label='SDSS Spectrum')


    plt.plot(wl_halpha, n_halpha_1, 'c-')
    plt.plot(wl_halpha, n_nii_6548_1, 'c-')
    plt.plot(wl_halpha, n_nii_6584_1, 'c-')
    plt.plot(wl_halpha, n_sii_6717_1, 'c-')
    plt.plot(wl_halpha, n_sii_6731_1, 'c-')

#    plt.plot(wl_halpha, n_halpha_2, 'g-')
#    plt.plot(wl_halpha, n_nii_6548_2, 'g-')
#    plt.plot(wl_halpha, n_nii_6584_2, 'g-')
#    plt.plot(wl_halpha, n_sii_6717_2, 'g-')
#    plt.plot(wl_halpha, n_sii_6731_2, 'g-')
#
#    plt.plot(wl_halpha, b_halpha_1, 'k-')

    ab = []
    for i in range(len(wl_halpha)):
        ab.append(m3(wl_halpha[i])+np.random.normal(loc=np.median(fl_halpha), scale=1.0))

#    plt.plot(wl_halpha, m3(wl_halpha), 'r-', label='MRES Spectrum')
    plt.plot(wl_halpha, ab, 'r-', label='MRES Spectrum')
    plt.xlim(6500, 6800)
    
    plt.xlabel('Wavelength ($\\rm \AA$)', fontsize='x-large')
    plt.ylabel('Flux Density ($\\rm \\times 10^{-17} \ erg \ cm^{-2} s^{-1} \AA^{-1}$)', fontsize='x-large')
    #plt.title('Emission Lines Fitting\nspec-' + plate[j] + '-' + mjd[j] + '-' + fiberid[j] + '.fits')
    #plt.savefig('figure3b/fig3b_'+ str(j) +'_spec-' + plate[j] + '-' + mjd[j] + '-' + fiberid[j])
    plt.savefig('fig3b_'+ str(j) +'_spec-' + plate[j] + '-' + mjd[j] + '-' + fiberid[j])
    plt.legend(loc='best')
    plt.show()
    
    
    print 'j =', j
#    plt.close('all')

    q+=1
#    if j % 50 == 0:
#	print 'Pausing for 5 seconds'
#	time.sleep(300.0)
