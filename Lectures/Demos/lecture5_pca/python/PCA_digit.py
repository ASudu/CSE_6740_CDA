############################################################################
# In this example we reduce the digit data(16X16) into 2 dimension 
# using PCA.
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
import show_image

# Loading the upsc digit dataset
matFile = sio.loadmat('usps_all.mat')

data = matFile['data']
pixelno = data.shape[0]
digitno = data.shape[1]
classno = data.shape[2]

# Displaying the digit 1(data(:,:,1)) and 0(data(:,:,10)) in the dataset
H = 16
W = 16
plt.figure(1)
digits_01 = np.concatenate((np.array(data[:,:,0]), np.array(data[:,:,9])), axis = 1).T
show_image.show_image_function(digits_01, H, W)
plt.title('digit 1 and 0')
#plt.figure(1)
#show_image.show_image_function(np.array(data[:,:,9]).T, H, W)
#plt.title('digit 2')

# Extracting the digits 1 and 0 and converting into double
x0 = np.concatenate((np.array(data[:, :, 0]), np.array(data[:, :, 9])), axis = 1) 
x = np.array(x0, dtype=np.float)
y = np.concatenate((np.ones((1,digitno)), 2 * np.ones((1,digitno))), axis = 1)

# number of data points to work with; 
m = x.shape[1]

# Normalize the data
Anew = x.T
stdA = np.std(Anew, axis = 0)
Anew = Anew.dot( np.diag(1 / stdA) )
Anew = Anew.T

# PCA
# Subtracting the mean of the dataset
mu = Anew.sum(axis=1) / m
mu  = mu.reshape((mu.shape[0],1))

xc = Anew - mu

# Finding the covariance 
C = xc.dot(xc.T) / m

# Finding top 2 pricipal component(eigen vector of the covariance)
k = 2
S, V = np.linalg.eig(C)
sortidx = S.argsort()[::-1] 
S = S[sortidx][0:k]
V = V[:,sortidx][:,0:k]
S = np.diag(S)

# python obtains the opposite eigenvector
V[:, 0] = -V[:, 0] 

# projecting the 16X16 data on the 2 eigenvectors
dim1 = V[:, 0].T.dot(xc) / np.sqrt(S[0,0]) 
dim2 = V[:, 1].T.dot(xc) / np.sqrt(S[1,1]) 

# Displaying the the data along the 2 PCs
plt.figure(2)
plt.hold(True)
plt.plot(dim1[(y == 1).squeeze()], dim2[(y == 1).squeeze()], 'r.')
plt.plot(dim1[(y == 2).squeeze()], dim2[(y == 2).squeeze()], 'b.')
plt.hold(False)
plt.show()