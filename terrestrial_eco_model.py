import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def lm(X,Xt,y,weights,n):
  W      = weights*np.eye(n)
  Sinv   = Xt@W@X
  L      = np.linalg.cholesky(Sinv)
  XtWy   = Xt@W@y
  temp   = np.linalg.solve(L, XtWy)
  B      = np.linalg.solve(np.transpose(L), temp)
  y_XB   = y-X@B
  sq_res = (y_XB)**2.
  sigma  = np.sqrt(np.transpose(y_XB)@W@y_XB/n)
  return([sq_res,sigma,B])

def design_matrix(time, harmonics, period,n):
  ones = np.ones(n)
  X = np.column_stack((ones,time))
  for i in range(harmonics):
    w = (2. * np.pi * time * (i+1) / period).reshape(n,1)
    X = np.concatenate((X,np.cos(w),np.sin(w)),axis=1)
  return(X)

def phenology_regression(time,y, nu, harmonics, period, maxit = 1000, epsilon = 1e-2, trace = True):
  n   = len(y)
  X   = design_matrix(time, harmonics, period,n) # initialize data
  Xt  = np.transpose(X)
  Q_t = 0          # complete data likelihood at iteration t
  Eq  = np.ones(n) # initialize latent variables
  #
  for it in range(maxit):
  # [ M-step]
    sq_res,sigma,B = lm(X, Xt, y, Eq,n)
    Q_tplus1       = np.sum(np.log(sigma) + (Eq*sq_res/(2*sigma**2))) # Compute the complete data log likelihood
  # [ E-step ]
    Eq = (nu+1)/(nu+((sq_res)/(sigma**2)))
    if trace: print("[", it, "] Q_t = ", Q_tplus1)
    if (it > 5) & ((Q_t-Q_tplus1) / abs(Q_t + .01) < epsilon):      # Check convergence
      break 
    Q_t = Q_tplus1
  return([B, sigma, Eq])


# Code for choosing fonts and colors
#
# from matplotlib import rcParams
# rcParams['font.family'] = 'Helvetica'
# rcParams['font.sans-serif'] = ['Helvetica Neue']
# from pylab import *
# cmap = cm.get_cmap('Purples')
# for i in range(cmap.N):
#     rgba = cmap(i)
#     # rgb2hex accepts rgb or rgba
#     print(matplotlib.colors.rgb2hex(rgba))

#Plot raw data
plt.scatter(time2,y,alpha=0.95,color='#3f007d')
plt.title('Raw data',fontsize=14)
plt.xlabel('Time',fontsize=14)
plt.ylabel('Spectral reflectance',fontsize=14)
plt.show()

# Fade outliers identified by the model
plt.scatter(time2,y,alpha=0.95,c=Eq,cmap='Purples')
plt.title('Data with weighted outliers',fontsize=14)
plt.xlabel('Time',fontsize=14)
plt.ylabel('Spectral reflectance',fontsize=14)
plt.show()

# Include the functional fit
time_pred = np.array([float(i/100.) for i in range(100*n_init)])
X_pred    = design_matrix(time_pred,harmonics,period,len(time_pred))
y_pred    = X_pred@B

plt.scatter(time2,y,alpha=0.95,c=Eq,cmap='Purples')
plt.title('Data with functional fit',fontsize=14)
plt.xlabel('Time',fontsize=14)
plt.ylabel('Spectral reflectance',fontsize=14)
plt.plot(time_pred,y_pred,color = 'black',linewidth=0.75)
plt.show()



