# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 14:05:40 2016

@author: Irham
"""


from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
from scipy.stats import linregress, spearmanr
import numpy as np
import pandas as pd

o = open('result/output_parameters.csv', 'w')
o.write('number\tm\tc\tr\tp\tstd\trho\tP\n')
o.close()

df2 = pd.read_csv('main_data.csv')
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

DP_oiii = DP.loc[DP.shift_OIII_4959_1.dropna().index]
DP_oiii = DP_oiii.loc[DP_oiii.FWHM_OIII_4959_1.dropna().index]
DP_oiii = DP_oiii.loc[DP_oiii.L_OIII_4959_1.dropna().index]
DP_oiii = DP_oiii.loc[DP_oiii.shift_OIII_4959_2.dropna().index]
DP_oiii = DP_oiii.loc[DP_oiii.FWHM_OIII_4959_2.dropna().index]
DP_oiii = DP_oiii.loc[DP_oiii.L_OIII_4959_2.dropna().index]

DP_oiii = DP.loc[DP.shift_nHbeta_1.dropna().index]
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

dc = [2666,2695,3388,23886,3488,3521,3887,4034,4592,4780,6900,7104,8743,8793,
10045,10073,11987,12009,12290,12292,12671,12692,13441,13533,13555,13595,13604,
13692,15990,15999,16856,16912,17111,17188,17711,17725,18132,18144,20175,20272]

dc_dp = []


ind = []
for i in DP.index:
    if i in DP_oiii.index\
    or i in DP_halpha.index:
        ind.append(i)

dp = DP.loc[ind]
DP_T1 = dp[(dp['double_peaked'] == True) & (dp.subtype=='AGN') & (dp.Type == 1)]
DP_T2 = dp[(dp['double_peaked'] == True) & (dp.subtype=='AGN') & (dp.Type == 2)]


for i in range(len(dc)):
#    if dc[i] in DP_halpha.index:
#        dc_dp.append(dc[i])
#    if dc[i] in DP_oiii.index:
#        dc_dp.append(dc[i])        
    if dc[i] in dp.index:
        dc_dp.append(dc[i])




###############################################################################

'''========================== SFR vs Redshift =========================='''

Lcont_bin2_T1 = [[], [], [], [], [], [], [], [], [], []]
z_bin_T1 = [[], [], [], [], [], [], [], [], [], []]
log_z_1_T1 = [[], [], [], [], [], [], [], [], [], []]
SFR_bin_T1 = [[], [], [], [], [], [], [], [], [], []]

Lcont_bin2_T2 = [[], [], [], [], [], [], [], [], [], []]
z_bin_T2 = [[], [], [], [], [], [], [], [], [], []]
log_z_1_T2 = [[], [], [], [], [], [], [], [], [], []]
SFR_bin_T2 = [[], [], [], [], [], [], [], [], [], []]


Lcont_bin_T1 = np.histogram(T1['L_5100'].dropna())
Lcont_T1 = T1['L_5100']
z_T1 = T1['z']

Lcont_bin_T2 = np.histogram(T2['L_5100'].dropna())
Lcont_T2 = T2['L_5100']
z_T2 = T2['z']




for i in range(len(Lcont_bin_T1[1]) - 1):
    
    for j in Lcont_T1.index:#range(len(Lcont)):
        if Lcont_T1[j] > Lcont_bin_T1[1][i] and Lcont_T1[j] <= Lcont_bin_T1[1][i+1]:
                        
            if str(T1.SFR[j]) != str(np.nan):

                SFR_bin_T1[i].append(T1['SFR'][j])
                log_z_1_T1[i].append(np.log10(z_T1[j]+1))
                
                z_bin_T1[i].append(z_T1[j])            
                Lcont_bin2_T1[i].append(Lcont_T1[j])

        else:
            continue
#---------------------------------------------------------------
    for j in Lcont_T2.index:#range(len(Lcont)):
        if Lcont_T2[j] > Lcont_bin_T2[1][i] and Lcont_T2[j] <= Lcont_bin_T2[1][i+1]:
                        
            if str(T2.SFR[j]) != str(np.nan):

                SFR_bin_T2[i].append(T2['SFR'][j])
                log_z_1_T2[i].append(np.log10(z_T2[j]+1))
                
                z_bin_T2[i].append(z_T2[j])            
                Lcont_bin2_T2[i].append(Lcont_T2[j])

        else:
            continue

z_bin_T1 = np.array(z_bin_T1)
Lcont_bin2_T1 = np.array(Lcont_bin2_T1)
SFR_bin_T1 = np.array(SFR_bin_T1)
log_z_1_T1 = np.array(log_z_1_T1)

z_bin_T2 = np.array(z_bin_T2)
Lcont_bin2_T2 = np.array(Lcont_bin2_T2)
SFR_bin_T2 = np.array(SFR_bin_T2)
log_z_1_T2 = np.array(log_z_1_T2)

        
bin_z = np.arange(0., 0.35 + 0.05/1.5, 0.05/1.5)#min(z), max(z), 0.07)#0.025*1.5)
bin_log_z_1 = np.arange(0.01, 0.13, 0.01)


def plot_average2(lab, col, number, x, y, group, labelx, labely):
    
    binned = []    
    mean = []
    stdx = []
    stdy = []
    count = []
    
    for i in range(len(group)):
        point = []
        axis = []
        for j in range(len(x)):
            if abs(x[j]-group[i]) <= abs((group[0]-group[1])/2.): 
                axis.append(x[j])
                point.append(y[j])
        
        if len(point) >= 5:
            if number >= 3000 and number <= 3006:
                mean.append(np.mean(point))
            else:
                mean.append(np.median(point))
                
            if number == 3000 or number == 3001:
                stdy.append(np.std(point)*1.)
            else:
                stdy.append(np.std(point)*0.25)
                
            stdx.append(np.std(axis)*1.)            
            count.append(len(point))
            binned.append(group[i])
    
    m, c, r, p, std = linregress(binned, mean)
    rho, P = spearmanr(binned, mean)
    o = open('result/output_parameters.csv', 'a')    
    o.write(str(number) + '\t' + str(m) + '\t' + str(c) + '\t' + str(r) + '\t' + str(p) + '\t' + str(std) + '\t' + str(rho) + '\t' + str(P) + '\n')
    o.close()
    plt.figure(number)
    plt.errorbar(binned, mean, yerr=stdy, xerr=0, fmt=str(col), label=str(lab))
    if number <= 3001:
        plt.plot(binned, mean, str(col)[0]+'-')
        

    plt.xlabel(labelx, fontsize='x-large')
    plt.ylabel(labely, fontsize='x-large')
    if number <= 3006:
        plt.xlim(min(binned)-0.15, max(binned)+0.05)
#    if number == 3006:
#        plt.xlim(-0.5, 10.)
#    if number == 3007:
#        plt.xlim(0, 20000)
        
    plt.legend(loc='lower right', fontsize='x-small')
    plt.savefig('figures/fig_%i' %number)
    
#plot_average2('L3', 'go', 3000, zO3_bin[2], (Lcont_bin2[2]), bin_z, 'Redshift $(z)$', '$\\log \ L_{5100}$')
#plot_average2('L4', 'bo', 3000, z_bin_T1[3], (Lcont_bin2_T1[3]), bin_z, 'Redshift $(z)$', '$\\log \ L_{5100}$')
plot_average2('L5', 'ko', 3000, z_bin_T1[4], (Lcont_bin2_T1[4]), bin_z, 'Redshift $(z)$', '$\\log \ L_{5100}$')
plot_average2('L6', 'co', 3000, z_bin_T1[5], (Lcont_bin2_T1[5]), bin_z, 'Redshift $(z)$', '$\\log \ L_{5100}$')
plot_average2('L7', 'yo', 3000, z_bin_T1[6], (Lcont_bin2_T1[6]), bin_z, 'Redshift $(z)$', '$\\log \ L_{5100}$')
plot_average2('L8', 'mo', 3000, z_bin_T1[7], (Lcont_bin2_T1[7]), bin_z, 'Redshift $(z)$', '$\\log \ L_{5100}$')
#plot_average2('L9', 'ro', 3000, zO3_bin[8], (Lcont_bin2[8]), bin_z, 'Redshift $(z)$', '$\\log \ L_{5100}$')


#plot_average2('L4', 'bo', 3005, log_z_1_T1[3], (SFR_bin_T1[3]), bin_log_z_1, '$\\log \ (z+1)$', '$\\log \ \\rm SFR \ M_{\\odot} \ yr^{-1}$')
plot_average2('L5', 'ko', 3005, log_z_1_T1[4], (SFR_bin_T1[4]), bin_log_z_1, '$\\log \ (z+1)$', '$\\log \ \\rm SFR \ M_{\\odot} \ yr^{-1}$')
plot_average2('L6', 'co', 3005, log_z_1_T1[5], (SFR_bin_T1[5]), bin_log_z_1, '$\\log \ (z+1)$', '$\\log \ \\rm SFR \ M_{\\odot} \ yr^{-1}$')
plot_average2('L7', 'yo', 3005, log_z_1_T1[6], (SFR_bin_T1[6]), bin_log_z_1, '$\\log \ (z+1)$', '$\\log \ \\rm SFR \ M_{\\odot} \ yr^{-1}$')
plot_average2('L8', 'mo', 3005, log_z_1_T1[7], (SFR_bin_T1[7]), bin_log_z_1, '$\\log \ (z+1)$', '$\\log \ \\rm SFR \ M_{\\odot} \ yr^{-1}$')

plt.figure(3005)
#plt.plot(bin_log_z_1, linear(bin_log_z_1, 2.02, 0.68), 'r-', label='$0.68 + 2.02 \ \\log(z+1)$')
#plt.plot(bin_log_z_1, linear(bin_log_z_1, 4.50, 0.55), 'k-', label='$0.55 + 4.50 \ \\log(z+1)$')
plt.xlim(min(bin_log_z_1)-0.01, max(bin_log_z_1)+0.01)
#plt.legend(loc='lower right', fontsize='x-small')
plt.savefig('figures/fig_3005')


#plot_average2('L3', 'go', 3000, zO3_bin[2], (Lcont_bin2[2]), bin_z, 'Redshift $(z)$', '$\\log \ L_{5100}$')
#plot_average2('L4', 'bo', 3001, z_bin_T2[3], (Lcont_bin2_T2[3]), bin_z, 'Redshift $(z)$', '$\\log \ L_{5100}$')
plot_average2('L5', 'ko', 3001, z_bin_T2[4], (Lcont_bin2_T2[4]), bin_z, 'Redshift $(z)$', '$\\log \ L_{5100}$')
plot_average2('L6', 'co', 3001, z_bin_T2[5], (Lcont_bin2_T2[5]), bin_z, 'Redshift $(z)$', '$\\log \ L_{5100}$')
plot_average2('L7', 'yo', 3001, z_bin_T2[6], (Lcont_bin2_T2[6]), bin_z, 'Redshift $(z)$', '$\\log \ L_{5100}$')
plot_average2('L8', 'mo', 3001, z_bin_T2[7], (Lcont_bin2_T2[7]), bin_z, 'Redshift $(z)$', '$\\log \ L_{5100}$')
#plot_average2('L9', 'ro', 3000, zO3_bin[8], (Lcont_bin2[8]), bin_z, 'Redshift $(z)$', '$\\log \ L_{5100}$')

#plot_average2('L4', 'bo', 3006, log_z_1_T2[3], (SFR_bin_T2[3]), bin_log_z_1, '$\\log \ (z+1)$', '$\\log \ \\rm SFR \ M_{\\odot} \ yr^{-1}$')
plot_average2('L5', 'ko', 3006, log_z_1_T2[4], (SFR_bin_T2[4]), bin_log_z_1, '$\\log \ (z+1)$', '$\\log \ \\rm SFR \ M_{\\odot} \ yr^{-1}$')
plot_average2('L6', 'co', 3006, log_z_1_T2[5], (SFR_bin_T2[5]), bin_log_z_1, '$\\log \ (z+1)$', '$\\log \ \\rm SFR \ M_{\\odot} \ yr^{-1}$')
plot_average2('L7', 'yo', 3006, log_z_1_T2[6], (SFR_bin_T2[6]), bin_log_z_1, '$\\log \ (z+1)$', '$\\log \ \\rm SFR \ M_{\\odot} \ yr^{-1}$')
plot_average2('L8', 'mo', 3006, log_z_1_T2[7], (SFR_bin_T2[7]), bin_log_z_1, '$\\log \ (z+1)$', '$\\log \ \\rm SFR \ M_{\\odot} \ yr^{-1}$')

plt.figure(3006)
#plt.plot(bin_log_z_1, linear(bin_log_z_1, 2.02, 0.68), 'r-', label='$0.68 + 2.02 \ \\log(z+1)$')
#plt.plot(bin_log_z_1, linear(bin_log_z_1, 4.50, 0.55), 'k-', label='$0.55 + 4.50 \ \\log(z+1)$')
plt.xlim(min(bin_log_z_1)-0.01, max(bin_log_z_1)+0.01)
#plt.legend(loc='lower right', fontsize='x-small')
plt.savefig('figures/fig_3006')
plt.close('all')

###############################################################################


def Ka03(x):#x = log(NII/Halpha), y = log(OIII/Hbeta)
    return 0.61/(x-0.05) + 1.3

def Ke01_a(x):#x = log(NII/Halpha), y = log(OIII/Hbeta)
    return 0.61/(x-0.47) + 1.19

def Ke01_b(x):#x = log(SII/Halpha), y = log(OIII/Hbeta)
    return 0.72/(x-0.32) + 1.3

def Ke06(x):#x = log(SII/Halpha), y = log(OIII/Hbeta)
    return 1.89*x + 0.76


plt.figure(3)
plt.plot(dp['LNII_nHalpha'], dp['LOIII_nHbeta'], 'k.', alpha=0.5)
plt.plot(np.linspace(-3, 0.2, 100), Ke01_a(np.linspace(-3, 0.2, 100)), 'r-', linewidth=3)
plt.plot(np.linspace(-3, -0.2, 100), Ka03(np.linspace(-3, -0.2, 100)), 'b-', linewidth=3)
plt.xlabel('$\\log \ \\rm [N \ II]/H\\alpha $', fontsize='x-large')
plt.ylabel('$\\log \ \\rm [O \ III]/H\\beta $', fontsize='x-large')
plt.ylim(-1.1, 2)
plt.xlim(-1.5, 0.8)
plt.title('Diagnostic Diagram', fontsize='x-large')
plt.savefig('figures/fig_3_bpt_nii_dp')

plt.figure(4)
plt.plot(dp['LSII_nHalpha'], dp['LOIII_nHbeta'], 'k.', alpha=0.5)
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

######################### FAILED #########################
'''
plt.figure(7)
#plt.plot(abs(dp.FWHM_nHalpha_1 - dp.FWHM_nHalpha_2)\
#    /dp[['FWHM_nHalpha_1','FWHM_nHalpha_2']].min(axis=1),\
#    (dp.L_nHalpha_1 - dp.L_nHalpha_2),\
#     'bo')
plt.plot(abs(dp.FWHM_OIII_5007_1 - dp.FWHM_OIII_5007_2)\
    /dp[['FWHM_OIII_5007_1','FWHM_OIII_5007_2']].min(axis=1),\
    (dp.L_OIII_5007_1 - dp.L_OIII_5007_2),\
     'bo')

plt.xlim(0.8, 5)
brak

dp_oiii = dp  [
            ( (abs(dp.FWHM_OIII_5007_1 - dp.FWHM_OIII_5007_2)
            /dp[['FWHM_OIII_5007_1','FWHM_OIII_5007_2']].min(axis=1)) > 0.8 )
            &
            ( (abs(dp.FWHM_OIII_5007_1 - dp.FWHM_OIII_5007_2)
            /dp[['FWHM_OIII_5007_1','FWHM_OIII_5007_2']].min(axis=1)) < 5. )    
            ]

dp_halpha = dp  [
            ( (abs(dp.FWHM_nHalpha_1 - dp.FWHM_nHalpha_2)
            /dp[['FWHM_nHalpha_1','FWHM_nHalpha_2']].min(axis=1)) > 0.8 )
            &
            ( (abs(dp.FWHM_nHalpha_1 - dp.FWHM_nHalpha_2)
            /dp[['FWHM_nHalpha_1','FWHM_nHalpha_2']].min(axis=1)) < 5. )
            ]
    
'''

'''
shift_r = pd.Series(np.zeros_like(DP_T1.shift_OIII_5007_1), index=DP_T1.index).replace(0, np.nan)
shift_b = pd.Series(np.zeros_like(DP_T1.shift_OIII_5007_1), index=DP_T1.index).replace(0, np.nan)

for i in DP_T1.shift_OIII_5007_1.index:
    if DP_T1.shift_OIII_5007_1.loc[i] >= DP_T1.shift_OIII_5007_2.loc[i]:   
        shift_r.loc[i] = DP_T1.shift_OIII_5007_1.loc[i]
        shift_b.loc[i] = DP_T1.shift_OIII_5007_2.loc[i]        
    else:
        shift_r.loc[i] = DP_T1.shift_OIII_5007_2.loc[i]
        shift_b.loc[i] = DP_T1.shift_OIII_5007_1.loc[i]
        
plt.figure(4)
plt.hist(shift_r.dropna(), bins=30, alpha=0.5)
plt.hist(shift_b.dropna(), bins=30, alpha=0.5)
'''


sel_oiii = np.loadtxt('notes/list_sel_test.csv', delimiter=',', dtype=int)
inde = []

for i in DP_oiii.index:
    if i in sel_oiii:
        inde.append(i)

dp = DP_oiii.loc[inde] #parameter ini cek lagi
shift_r = pd.Series(np.zeros_like(dp.shift_OIII_5007_1), index=dp.index).replace(0, np.nan)
shift_b = pd.Series(np.zeros_like(dp.shift_OIII_5007_1), index=dp.index).replace(0, np.nan)

l_r = pd.Series(np.zeros_like(dp.L_OIII_5007_1), index=dp.index).replace(0, np.nan)
l_b = pd.Series(np.zeros_like(dp.L_OIII_5007_1), index=dp.index).replace(0, np.nan)

fwhm_r = pd.Series(np.zeros_like(dp.FWHM_OIII_5007_1), index=dp.index).replace(0, np.nan)
fwhm_b = pd.Series(np.zeros_like(dp.FWHM_OIII_5007_1), index=dp.index).replace(0, np.nan)

delta_rb = pd.Series(np.zeros_like(dp.shift_OIII_5007_1), index=dp.index).replace(0, np.nan)
fwhm_min = pd.Series(np.zeros_like(dp.shift_OIII_5007_1), index=dp.index).replace(0, np.nan)



dp_halpha = DP_halpha

shift_r_halpha = pd.Series(np.zeros_like(dp_halpha.shift_nHalpha_1), index=dp_halpha.index).replace(0, np.nan)
shift_b_halpha = pd.Series(np.zeros_like(dp_halpha.shift_nHalpha_1), index=dp_halpha.index).replace(0, np.nan)

l_r_halpha = pd.Series(np.zeros_like(dp_halpha.L_nHalpha_1), index=dp_halpha.index).replace(0, np.nan)
l_b_halpha = pd.Series(np.zeros_like(dp_halpha.L_nHalpha_1), index=dp_halpha.index).replace(0, np.nan)

fwhm_r_halpha = pd.Series(np.zeros_like(dp_halpha.FWHM_nHalpha_1), index=dp_halpha.index).replace(0, np.nan)
fwhm_b_halpha = pd.Series(np.zeros_like(dp_halpha.FWHM_nHalpha_1), index=dp_halpha.index).replace(0, np.nan)

delta_rb_halpha = pd.Series(np.zeros_like(dp_halpha.shift_nHalpha_1), index=dp_halpha.index).replace(0, np.nan)
fwhm_min_halpha = pd.Series(np.zeros_like(dp_halpha.shift_nHalpha_1), index=dp_halpha.index).replace(0, np.nan)

#delta_rb = shift_r-shift_b
#fwhm_min = dp[['FWHM_nHalpha_1','FWHM_nHalpha_2']].min(axis=1)


for i in dp.shift_OIII_5007_1.index:
    if str(dp.shift_OIII_5007_1.loc[i]) == str(np.nan)\
    or str(dp.shift_OIII_5007_2.loc[i]) == str(np.nan):
        continue
    
    if dp.shift_OIII_5007_1.loc[i] >= dp.shift_OIII_5007_2.loc[i]:   
        shift_r.loc[i] = dp.shift_OIII_5007_1.loc[i]
        shift_b.loc[i] = dp.shift_OIII_5007_2.loc[i]

        l_r.loc[i] = dp.L_OIII_5007_1.loc[i]
        l_b.loc[i] = dp.L_OIII_5007_2.loc[i]

        fwhm_r.loc[i] = dp.FWHM_OIII_5007_1.loc[i]
        fwhm_b.loc[i] = dp.FWHM_OIII_5007_2.loc[i]
        
        delta_rb.loc[i] = shift_r.loc[i] - shift_b.loc[i]
        fwhm_min.loc[i] = np.min([fwhm_r.loc[i], fwhm_b.loc[i]])
        
    else:
        shift_r.loc[i] = dp.shift_OIII_5007_2.loc[i]
        shift_b.loc[i] = dp.shift_OIII_5007_1.loc[i]

        l_r.loc[i] = dp.L_OIII_5007_2.loc[i]
        l_b.loc[i] = dp.L_OIII_5007_1.loc[i]

        fwhm_r.loc[i] = dp.FWHM_OIII_5007_2.loc[i]
        fwhm_b.loc[i] = dp.FWHM_OIII_5007_1.loc[i]

        delta_rb.loc[i] = shift_r.loc[i] - shift_b.loc[i]
        fwhm_min.loc[i] = np.min([fwhm_r.loc[i], fwhm_b.loc[i]])

        
        
for i in dp_halpha.shift_nHalpha_1.index:
    if str(dp_halpha.shift_nHalpha_1.loc[i]) == str(np.nan)\
    or str(dp_halpha.shift_nHalpha_2.loc[i]) == str(np.nan):
        continue
    
    if dp_halpha.shift_nHalpha_1.loc[i] >= dp_halpha.shift_nHalpha_2.loc[i]:   
        shift_r_halpha.loc[i] = dp_halpha.shift_nHalpha_1.loc[i]
        shift_b_halpha.loc[i] = dp_halpha.shift_nHalpha_2.loc[i]

        l_r_halpha.loc[i] = dp_halpha.L_NII_6584_1.loc[i]
        l_b_halpha.loc[i] = dp_halpha.L_NII_6584_2.loc[i]

        fwhm_r_halpha.loc[i] = dp_halpha.FWHM_nHalpha_1.loc[i]
        fwhm_b_halpha.loc[i] = dp_halpha.FWHM_nHalpha_2.loc[i]
        
        delta_rb_halpha.loc[i] = shift_r_halpha.loc[i] - shift_b_halpha.loc[i]
        fwhm_min_halpha.loc[i] = np.min([fwhm_r_halpha.loc[i], fwhm_b_halpha.loc[i]])
        
    else:
        shift_r_halpha.loc[i] = dp_halpha.shift_nHalpha_2.loc[i]
        shift_b_halpha.loc[i] = dp_halpha.shift_nHalpha_1.loc[i]

        l_r_halpha.loc[i] = dp_halpha.L_NII_6584_2.loc[i]
        l_b_halpha.loc[i] = dp_halpha.L_NII_6584_1.loc[i]

        fwhm_r_halpha.loc[i] = dp_halpha.FWHM_nHalpha_2.loc[i]
        fwhm_b_halpha.loc[i] = dp_halpha.FWHM_nHalpha_1.loc[i]

        delta_rb_halpha.loc[i] = shift_r_halpha.loc[i] - shift_b_halpha.loc[i]
        fwhm_min_halpha.loc[i] = np.min([fwhm_r_halpha.loc[i], fwhm_b_halpha.loc[i]])


#coba buat plot broooo

log_v = np.median(np.log10(abs(shift_r/shift_b)))
plt.figure(11)
plt.hist(shift_b, edgecolor='b', linewidth=2, facecolor='none', bins=20, label='$\\langle \\Delta V_b \\rangle = %.2f$' %(np.median(shift_b.dropna())))
plt.hist(shift_r, edgecolor='r', linewidth=2, facecolor='none', bins=20, label='$\\langle \\Delta V_r \\rangle = %.2f$' %(np.median(shift_r.dropna())))
plt.xlabel('$\\Delta V_{r, b} \ (\\rm km \ s^{-1})$', fontsize='x-large')
plt.ylabel('Number of AGN', fontsize='x-large')
plt.legend(loc='best')
plt.savefig('figures/fig_11')
plt.close('all')


ge = (delta_rb/fwhm_min > 0.8) & (delta_rb/fwhm_min < 3.4)
plt.plot(delta_rb/fwhm_min[ge], l_r-l_b[ge], 'bo')
plt.xlim(0, 4)
plt.close('all')

'''
plt.hist(np.log10(abs(shift_r/shift_b)))
plt.hist((l_r-l_b).dropna())
plt.xlim(-1.5, 1.5)
plt.ylim(-1.5, 1.5)
plt.plot(np.log10(abs(shift_r/shift_b)), (l_r-l_b), 'bo')
plt.xlabel('$\\rm \\log \ (V_b/V_r)$', fontsize='x-large')
plt.ylabel('$\\rm \\log \ (F_b/F_r)$', fontsize='x-large')
plt.title('$\\rm [O III] \ \lambda5007$')
plt.legend(loc='best')
plt.close('all')
'''

plt.plot(np.log10(abs(shift_r_halpha/shift_b_halpha)), (l_r_halpha-l_b_halpha), 'bo')
plt.xlabel('$\\rm \\log \ (V_b/V_r)$', fontsize='x-large')
plt.ylabel('$\\rm \\log \ (F_b/F_r)$', fontsize='x-large')
plt.title('$\\rm H\\alpha$')
plt.legend(loc='best')
#plt.xlim(-1.5, 1.5)
#plt.ylim(-1.5, 1.5)