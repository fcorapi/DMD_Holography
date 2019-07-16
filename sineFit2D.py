#sineFit2D
#This program fits a 2-dimensional sine function to an image. Used for determing the phase parameters
#for the DMD interference patterns in order to obtain a calibration phase map.
#Modified the 2D Gaussian Fit code from Scipy Cookbook (https://scipy-cookbook.readthedocs.io/items/FittingData.html)
#Frank Corapi (fcorapi@uwaterloo.ca)
#July 16th, 2019

import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.image as img
from scipy import optimize

#Define the fitting function, a 2-dimensional sine function
def sine2D(A, B, x_0, y_0, phi):
    A = float(A)
    B = float(B)
    x_0 = float(x_0)
    y_0 = float(y_0)
    phi = float(phi)
    return lambda x,y: A**2 + B**2 + 2*A*B*np.cos(x_0*x + y_0*y + phi)

#Determines the indices of the local minima within a list
def local_min(list):
    minIndices = [i for i,y in enumerate(list)
                  if ((i == 0) or (list[i-1] >= y))
                  and ((i == len(list)-1) or (y < list[i+1]))]
    return minIndices

#This function is used to make an accurate initial guess for the fitting paramters for the 2D sine function
def guessParams(data):
    amplitude = data.max()

    A = np.sqrt(amplitude)/2
    B = A

    row = data[int(np.shape(data)[0]/2), :]
    col = data[:, int(np.shape(data)[1]/2)]
    minIndicesX = local_min(col)
    minIndicesY = local_min(row)

    periodGuessX = 0
    nx = 0
    for loop in range(0,len(minIndicesX)-1):
        guessX = minIndicesX[loop+1] - minIndicesX[loop]
        periodGuessX = periodGuessX + guessX
        nx = nx+1
    if nx == 0:
        x_0 = 0
    else:
        periodGuessX = periodGuessX/nx
        x_0 = 2*np.pi/periodGuessX

    periodGuessY = 0
    ny = 0
    for loop in range(0, len(minIndicesY) - 1):
        guessY = minIndicesY[loop + 1] - minIndicesY[loop]
        periodGuessY = periodGuessY + guessY
        ny = ny + 1

    if ny == 0:
        y_0 = 0
    else:
        periodGuessY = periodGuessY / ny
        y_0 = 2*np.pi/periodGuessY

    # print A,B,2*np.pi/x_0, 2*np.pi/y_0
    phi = 0

    return A, B, x_0, -y_0, phi

#Using the initial guess parameters from guessParams, the optimal parameters are determined using the following function
#A least Squares method is used.
def fitSine2D(data):
    params = guessParams(data)
    paramBounds = ([0,0,0,-np.inf,0],[np.inf,np.inf,np.inf,np.inf,2*np.pi])
    errorfunction = lambda p: np.ravel(sine2D(*p)(*np.indices(data.shape)) -
                                 data)
    p = optimize.least_squares(errorfunction, params, bounds = paramBounds)
    # print success
    return p.x

#Image Directories
dir = 'C:\Users\Franky\Desktop\UofT Summer 2019\CalibrationImages\\'
cropDir = 'C:\Users\Franky\Desktop\UofT Summer 2019\CalibrationImages\Cropped\\'
filename = 'CI_X6Y3'
filenamePrefix = 'CI_'
ext1 = '.bmp'
#Select the cropping region
cropCoords = (600,680,690,840)

image = Image.open(dir+filename+ext1)
croppedImage = image.crop(cropCoords)
croppedImage.show()
croppedData = np.array(croppedImage)[:,:,0]

# croppedImage.show()
plt.imshow(croppedData)
# plt.show()
# data = np.array(image)[:,:,0]
# plt.imshow(data)
# plt.show()

params = fitSine2D(croppedData)
print params
print 2*np.pi/params[2], 2*np.pi/params[3]
fit = sine2D(*params)

plt.contour(fit(*np.indices(croppedData.shape)), cmap=plt.cm.copper)
plt.show()
