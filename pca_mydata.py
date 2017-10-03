# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 13:41:25 2016

@author: irham
"""





import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

import numpy as np
from sklearn.preprocessing import normalize
import pandas as pd

#new_data = pd.read_csv('for_pca_mydata_halpha.csv', delimiter=',').dropna()
#
#
##new_data = pd.read_csv('for_pca_mydata_oiii.csv', delimiter=',').dropna()
##new_data.to_csv('for_pca_mydata_oiii_clean.csv')
##new_data.to_csv('for_pca_clean_classifier.csv')
#new_data.to_csv('test.csv')



data = np.loadtxt('for_pca_mydata_halpha_clean.csv', skiprows=1, delimiter = ',')
#data = np.loadtxt('for_pca_mydata_oiii_clean.csv', skiprows=1, delimiter = ',')
#data = np.loadtxt('result_for_pca_b.csv', skiprows=1, delimiter = ',')
#data = normalize(data, norm='l2')



def norm(x):
    for j in range(len(x[0])):
        
        x[:, j] = (x[:, j] - np.mean(x[:, j])) / (np.std(x[:, j])*np.sqrt(len(x)))
    
    return x


data = norm(data)

pca2 = PCA(n_components=5)
data_r = pca2.fit_transform(data)

print('explained variance ratio (first two components): \n%s'
      % str(pca2.explained_variance_ratio_))

print '\n========================================\n\n'

print '\nVariance (Eigenvalue)'
print pca2.explained_variance_

print '\nVariance Ratio (Proportion)'
print pca2.explained_variance_ratio_

#print '\nMean'
#print pca2.mean_


print '\nv_disp, shift_b_oiii, shift_r_oiii, fwhm_b_oiii, fwhm_r_oiii, l_b_oiii, l_r_oiii, l_r_oiii/l_b_oiii	, shift_r_oiii/shift_b_oiii'
print '\nComponents'
print pca2.components_

print '\nCumulative Variance (Cumulative)'
print pca2.explained_variance_ratio_.cumsum()


o = open('hasil.csv', 'w')

o.write('---\tEV1\tEV2\tEV3\tEV4\tEV5\n')
o.write('Eigenvalue\t')
for i in pca2.explained_variance_:
    o.write(str(i) + '\t')

o.write('\n')

o.write('Proportion\t')
for i in pca2.explained_variance_ratio_:
    o.write(str(i) + '\t')

o.write('\n')

o.write('Cumulative\t')
for i in pca2.explained_variance_ratio_.cumsum():
    o.write(str(i) + '\t')

o.write('\n')


name = ['v_disp', 'shift_b_oiii', 'shift_r_oiii', 'fwhm_b_oiii', 'fwhm_r_oiii', 'l_b_oiii', 'l_r_oiii', 'l_r_oiii/l_b_oiii', 'shift_r_oiii/shift_b_oiii']

j = 0

for i in range(len(name)):
    o.write(str(name[i]) + '\t')
    for j in range(5):
        o.write(str(pca2.components_[j][i]) + '\t')
    o.write('\n')

o.close()

#1 halpha, 2 oiii
index_dn = np.loadtxt('for_pca_classifier_index.csv', delimiter=',', skiprows=1, dtype=int)
masked_dn = np.loadtxt('masked.csv', delimiter=',', skiprows=1, dtype=int)
#index_of = np.loadtxt('for_pca_classifier_index_2.csv', delimiter=',', skiprows=1, dtype=int)

index_of = np.loadtxt('outflow_halpha_index.csv', delimiter=',', skiprows=1, dtype=int)

plt.figure(1)
plt.scatter(data_r[:, 0], data_r[:, 1])
plt.scatter(data_r[:, 0][index_dn], data_r[:, 1][index_dn], color='red', label = 'Binary AGN')
plt.scatter(data_r[:, 0][masked_dn], data_r[:, 1][masked_dn], color='green', label = 'Binary?')
plt.scatter(data_r[:, 0][index_of], data_r[:, 1][index_of], color='yellow', label = 'Outflow')
#plt.scatter(data_r[:, 0][index_of], data_r[:, 1][index_of], color='green')
plt.xlabel('Component 1', fontsize='x-large')
plt.ylabel('Component 2', fontsize='x-large')
plt.legend(loc='best')
plt.savefig('figures/fig_1')

plt.figure(2)
plt.scatter(data_r[:, 0], data_r[:, 2])
#plt.scatter(data_r[:, 0][np.arange(0, 1495, 3)], data_r[:, 2][np.arange(0, 1495, 3)])
plt.scatter(data_r[:, 0][index_dn], data_r[:, 2][index_dn], color='red', label = 'Binary AGN')
plt.scatter(data_r[:, 0][masked_dn], data_r[:, 2][masked_dn], color='green', label = 'Binary?')
plt.scatter(data_r[:, 0][index_of], data_r[:, 2][index_of], color='yellow', label = 'Outflow')
#plt.scatter(data_r[:, 0][index_of], data_r[:, 2][index_of], color='green')
plt.xlabel('Component 1', fontsize='x-large')
plt.ylabel('Component 3', fontsize='x-large')
plt.legend(loc='best')
plt.savefig('figures/fig_2')

plt.figure(3)
plt.scatter(data_r[:, 1], data_r[:, 2])
plt.scatter(data_r[:, 1][index_dn], data_r[:, 2][index_dn], color='red', label = 'Binary AGN')
plt.scatter(data_r[:, 1][masked_dn], data_r[:, 2][masked_dn], color='green', label = 'Binary?')
plt.scatter(data_r[:, 1][index_of], data_r[:, 2][index_of], color='yellow', label = 'Outflow')
#plt.scatter(data_r[:, 1][index_of], data_r[:, 2][index_of], color='green')
plt.xlabel('Component 2', fontsize='x-large')
plt.ylabel('Component 3', fontsize='x-large')
plt.legend(loc='best')
plt.savefig('figures/fig_3')


plt.close('all')

'''
Jadi gini. Hasil PCA untuk data Ge menunjukkan bahwa PCA belum mampu memisahkan
binary AGN dan outflow AGN dengan baik. Hal yang menarik adalah kebanyakan
binary berada di bagian tengah dari distribusi, sementara yang outflow cenderung
menjadi outlier
'''











#plt.close('all')

#from mpl_toolkits.mplot3d import Axes3D
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
##ax.scatter(data_r[:, 0], data_r[:, 1], data_r[:, 2])
#ax.scatter(data_r[:, 0][index_dn], data_r[:, 1][index_dn], data_r[:, 2][index_dn], c='red')
#ax.scatter(data_r[:, 0][index_of], data_r[:, 1][index_of], data_r[:, 2][index_of], c='green')
#plt.close('all')

#'''========== LLE =========='''
#
#from sklearn import manifold, neighbors
#from astroML.plotting.tools import discretize_cmap
#from astroML.decorators import pickle_results
#
##------------------------------------------------------------
## Set up color-map properties
#clim = (1.5, 6.5)
#cmap = discretize_cmap(plt.cm.jet, 5)
#cdict = ['DN, OF, U']
#cticks = [2, 3, 4, 5, 6]
#formatter = plt.FuncFormatter(lambda t, *args: cdict[int(np.round(t))])
#
#
##------------------------------------------------------------
## Compute the LLE projection; save the results
#
#spec = data
#
##@pickle_results("spec_LLE.pkl")
#def compute_spec_LLE(n_neighbors=10, out_dim=3):
#    # Compute the LLE projection
#    LLE = manifold.LocallyLinearEmbedding(n_neighbors, out_dim,
#                                          method='modified',
#                                          eigen_solver='dense')
#    Y_LLE = LLE.fit_transform(spec)
#    print " - finished LLE projection"
#
#    # remove outliers for the plot
#    BT = neighbors.BallTree(Y_LLE)
#    dist, ind = BT.query(Y_LLE, n_neighbors)
#    dist_to_n = dist[:, -1]
#    dist_to_n -= dist_to_n.mean()
#    std = np.std(dist_to_n)
#    flag = (dist_to_n > 0.25 * std)
#    print " - removing %i outliers for plot" % flag.sum()
#
#    return flag, Y_LLE#[~flag]#, color[~flag]
#
##coeffs_LLE, c_LLE = compute_spec_LLE(10, 3)
#fla, coeffs_LLE = compute_spec_LLE(10, 3)
#
#
#plt.figure(11)
#plt.scatter(coeffs_LLE[:, 0], coeffs_LLE[:, 1])
#plt.scatter(coeffs_LLE[:, 0][index_dn], coeffs_LLE[:, 1][index_dn], color='red')
#plt.scatter(coeffs_LLE[:, 0][index_of], coeffs_LLE[:, 1][index_of], color='green')
#plt.xlabel('Component 1')
#plt.ylabel('Component 2')