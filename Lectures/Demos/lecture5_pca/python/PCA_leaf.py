############################################################################
# In this example we visualize the leaf dataset which has 14 dimensions by
# transforming the data using PCA to project it along its 2 dominant basis.
# This method is often used for visualizing high dimentional data.
# leaf dataset
# ------------------------------
# The codes are based on Python2.7. 
# Please install numpy, scipy, matplotlib packages before using.
# Thank you for your suggestions!
#
# @version 1.0
# @author CSE6740/CS7641 TA
############################################################################
import scipy.io as sio 
import numpy as np
import matplotlib.pyplot as plt

data = sio.loadmat('leaf.mat')['M']

# dataset description
# The provided data comprises the following shape (attributes 3 to 9) and texture (attributes 10
# to 16) features:
# 1. Class (Species)
# 2. Specimen Number
# 3. Eccentricity
# 4. Aspect Ratio
# 5. Elongation
# 6. Solidity
# 7. Stochastic Convexity
# 8. Isoperimetric Factor
# 9. Maximal Indentation Depth
# 10. Lobedness
# 11. Average Intensity
# 12. Average Contrast
# 13. Smoothness
# 14. Third moment
# 15. Uniformity
# 16. Entropy


# extract attributes from the raw data
Anew = data[:,2:16];

m = Anew.shape[0]
n = Anew.shape[1]

# create indicator matrix; 
Inew = data[:,0]; 
# normalize data; 
stdA = np.std(Anew, axis = 0)
Anew = Anew.dot( np.diag(1 / stdA) )
Anew = Anew.T

# PCA
mu = Anew.sum(axis=1) / m
mu  = mu.reshape((mu.shape[0],1))

xc = Anew - mu

C = xc.dot(xc.T) / m

k = 2
S, V = np.linalg.eig(C)
sortidx = S.argsort()[::-1] 
S = S[sortidx][0:k]
V = V[:,sortidx][:,0:k]
S = np.diag(S)

# python obtains the opposite eigenvector
V[:, 0] = -V[:, 0] 
V[:, 1] = -V[:, 1] 

dim1 = V[:, 0].T.dot(xc) / np.sqrt(S[0,0]) 
dim2 = V[:, 1].T.dot(xc) / np.sqrt(S[1,1])

color_string = 'bgrmck'
marker_string = '.+*o'

plt.figure(1)
plt.hold(True)
for i in range(1, np.max(Inew.astype(int) + 1)):
    plt.plot(dim1[Inew == i], dim2[Inew == i],\
              color_string[i % 5] + marker_string[i % 4] )
plt.hold(False)
print V[:,0:2]

plt.show()