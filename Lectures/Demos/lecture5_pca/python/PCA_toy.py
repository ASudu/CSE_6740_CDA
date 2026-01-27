############################################################################
# In this example we first generate a 2 dimensional dataset which is
# distributed along an axis which is 45 degrees to the x-axis. We then
# apply PCA to show that the axis it generates is alinged along this axis.
# ------------------------------
# The codes are based on Python2.7. 
# Please install numpy, scipy, matplotlib packages before using.
# Thank you for your suggestions!
#
# @version 1.0
# @author CSE6740/CS7641 TA
############################################################################
import numpy as np
import matplotlib.pyplot as plt

G = np.diag(np.array([2, 0.5]))
R = np.array([[np.cos(np.pi/4), -np.sin(np.pi/4)],[np.sin(np.pi/4), np.cos(np.pi/4)]])
m = 1000
# Multiplying with G results in scaling of x-axis by 2 and y-axis by 0.5
# Multiplying with R results in rotation of the data by pi/4 rads
x = R.dot(G).dot(np.random.randn(2, m)) + np.tile(np.array([0.5, -0.5]),(m,1)).T
x = np.concatenate((x, R.dot(G).dot(np.random.randn(2, m)) + np.tile(np.array([-0.5, 0.5]),(m,1)).T), axis = 1)
y = np.concatenate((np.ones((1, m)), 2 * np.ones((1, m))), axis = 1)

iscolor = 1

# Plotting the orginal random data generated
plt.figure(1)
plt.hold(True)
if iscolor == 1:
    plt.plot(x[0,(y == 1).squeeze()], x[1,(y == 1).squeeze()], 'r.')
    plt.plot(x[0,(y == 2).squeeze()], x[1,(y == 2).squeeze()], 'b.')
else:
    plt.plot(x[0,:], x[1,:], 'b.')
plt.hold(False)

plt.figure(2)
plt.subplot(121)
plt.hist(x[0,:], bins=100) 
plt.subplot(122)
plt.hist(x[1,:], bins=100) 

plt.show(block = False)

# number of data points to work with; 
m = x.shape[1]

#raw_input('press key to continue ...')

# Normalize the the data
Anew = x.T 
stdA = np.std(Anew, axis = 0)
Anew = Anew.dot( np.diag(1 / stdA) )
Anew = Anew.T

# PCA
# Subtracting by the mean
mu = Anew.sum(axis=1) / m
mu  = mu.reshape((mu.shape[0],1))

xc = Anew - mu

# Finding Covariance matrix
C = xc.dot(xc.T) / m

# Finding eigenvectors of the Covarance matrix
k = 2
S, V = np.linalg.eig(C)
sortidx = S.argsort()[::-1] 
S = S[sortidx][0:k]
V = V[:,sortidx][:,0:k]
S = np.diag(S)

# Pojecting the data along the principal component
dim1 = V[:, 0].T.dot(xc) / np.sqrt(S[0,0]) 
dim2 = V[:, 1].T.dot(xc) / np.sqrt(S[1,1])

# Plotting the data projected along the principal component
plt.figure(3)
plt.hold(True)
if iscolor == 1:
    plt.plot(dim1[(y == 1).squeeze()], dim2[(y == 1).squeeze()], 'r.')
    plt.plot(dim1[(y == 2).squeeze()], dim2[(y == 2).squeeze()], 'b.')
else:
    plt.plot(dim1, dim2, 'b.')
plt.hold(False)

plt.figure(4)
plt.subplot(121)
plt.hist(dim1, bins=100) 
plt.subplot(122)
plt.hist(dim2, bins=100) 
plt.show()