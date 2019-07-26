#sineFit2D
#This program fits a 2-dimensional sine function to an image. Used for determing the phase parameters
#for the DMD interference patterns in order to obtain a calibration phase map.
#Modified the 2D Gaussian Fit code from Scipy Cookbook (https://scipy-cookbook.readthedocs.io/items/FittingData.html)
#Frank Corapi (fcorapi@uwaterloo.ca)
#July 17th, 2019

#Import Directories
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.image as img
from scipy import optimize


#****************************FUNCTIONS*****************************************************************8
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
def guessParams(data,xPos, yPos):
    amplitude = data.max()

    A = np.sqrt(amplitude)/2
    B = A

    #Optical System Parameters
    wavelength = 0.760 #microns
    dPatch = 0.0616 #cm
    cameraPixelSize = 5.2 #microns/pixel
    f1 = 10.0 #cm
    f2 = 4.0 #cm
    f3 = 30.0 #cm

    #Position of the Center patch (X=1 Y=1 is top left corner of DMD)
    xCen = 8
    yCen = 5

    deltaX = xCen - xPos
    if deltaX == 0:
        deltaX = 1e-3

    deltaY = yCen - yPos
    if deltaY == 0:
        deltaY = 1e-3

    if xPos < xCen and yPos < yCen:
        theta = np.arctan(deltaX/deltaY)
    elif xPos < xCen and yPos == yCen:
        theta = np.arctan(deltaX / deltaY)
    elif xPos < xCen and yPos > yCen:
        theta = np.arctan(deltaX / deltaY)
    elif xPos == xCen and yPos > yCen:
        theta = np.arctan(deltaX / deltaY) + np.pi
    elif xPos > xCen and yPos > yCen:
        theta = np.arctan(deltaX / deltaY) + np.pi
    elif xPos > xCen and yPos == yCen:
        theta = np.arctan(deltaX / deltaY)
    elif xPos > xCen and yPos < yCen:
        theta = np.arctan(deltaX / deltaY) + np.pi
    else:
        theta = np.arctan(deltaX/deltaY)

    Xi = np.arctan((f2/f1)*np.sqrt(deltaX**2 + deltaY**2)*dPatch/f3)
    print deltaX, deltaY, theta,Xi
    periodGuessX = (wavelength/cameraPixelSize)/(np.sin(Xi)*np.cos(theta))
    periodGuessY = (wavelength/cameraPixelSize)/(np.sin(Xi)*np.sin(theta))

    x_0 = 2*np.pi/periodGuessX
    y_0 = 2*np.pi/periodGuessY

    phi = 0

    return A, B, x_0, y_0, phi

#Using the initial guess parameters from guessParams, the optimal parameters are determined using the following function
#A least Squares method is used.
def fitSine2D(data, xPos, yPos):
    params = guessParams(data, xPos, yPos)
    paramBounds = ([0,0,-np.inf,-np.inf,0],[np.inf,np.inf,np.inf,np.inf,2*np.pi])
    errorfunction = lambda p: np.ravel(sine2D(*p)(*np.indices(data.shape)) -
                                 data)
    p = optimize.least_squares(errorfunction, params, bounds = paramBounds)
    # print success
    return p.x


#****************END OF FUNCTIONS**************************************************************************************

#Image Directories
dir = 'C:\Users\Franky\Desktop\UofT Summer 2019\CalibrationImages3 (July 25)\\'
cropDir = 'C:\Users\Franky\Desktop\UofT Summer 2019\CalibrationImages3 (July 25)\Cropped\\'
fitDir = 'C:\Users\Franky\Desktop\UofT Summer 2019\CalibrationImages3 (July 25)\Fitted\\'
filename = 'CI3_X4Y7'
filenamePrefix = 'CI3_'
ext1 = '.bmp'
ext2 = '.png'
#Select the cropping region
cropCoords = (590,420,770,580)
xMax = 16
yMax = 10

for xLoop in range(1,xMax+1):
    for yLoop in range(1,yMax+1):
        if xLoop == 8 and yLoop == 5:
            continue
        image = Image.open(dir + filenamePrefix + 'X' + str(xLoop) + 'Y' + str(yLoop) + ext1)
        # data = np.array(image)[:, :, 0]
        # plt.imshow(data)
        # plt.show()

        croppedImage = image.crop(cropCoords)
        croppedImage.save(cropDir + filenamePrefix + 'X' + str(xLoop) + 'Y' + str(yLoop) + '_Cropped'+ ext1)
        croppedData = np.array(croppedImage)[:, :, 0]

        params = fitSine2D(croppedData, xLoop, yLoop)
        # print params
        # print 2 * np.pi / params[2], 2 * np.pi / params[3]
        fit = sine2D(*params)
        plt.figure()
        plt.imshow(croppedData)
        plt.contour(fit(*np.indices(croppedData.shape)), cmap=plt.cm.copper)
        plt.colorbar()
        plt.savefig(fitDir + filenamePrefix + 'X' + str(xLoop) + 'Y' + str(yLoop) + '_Fit' + ext2)







#*************************OLD ATTEMPT**************************
# row = data[int(np.shape(data)[0]/2), :]
# col = data[:, int(np.shape(data)[1]/2)]
# minIndicesX = local_min(col)
# minIndicesY = local_min(row)
#
# periodGuessX = 0
# nx = 0
# for loop in range(0,len(minIndicesX)-1):
#     guessX = minIndicesX[loop+1] - minIndicesX[loop]
#     periodGuessX = periodGuessX + guessX
#     nx = nx+1
# if nx == 0:
#     x_0 = 0
# else:
#     periodGuessX = periodGuessX/nx
#     periodGuessX = -44.95
#     x_0 = 2*np.pi/periodGuessX
#
# periodGuessY = 0
# ny = 0
# for loop in range(0, len(minIndicesY) - 1):
#     guessY = minIndicesY[loop + 1] - minIndicesY[loop]
#     periodGuessY = periodGuessY + guessY
#     ny = ny + 1
#
# if ny == 0:
#     y_0 = 0
# else:
#     periodGuessY = periodGuessY / ny
#     periodGuessY = 10000
#     y_0 = 2*np.pi/periodGuessY

# print A,B,2*np.pi/x_0, 2*np.pi/y_0


# ft = np.fft.fft2(data)
# ftShift = np.fft.fftshift(ft)
# magnitude = 200 * np.log(np.abs(ftShift))
# # plt.imshow(magnitude)
# # plt.show()
# height = np.shape(magnitude)[0]
# width = np.shape(magnitude)[1]
# # plt.imshow(magnitude[0:height/2-3,0:width/2-3])
# # plt.show()
# # plt.imshow(magnitude[0:height/2-3,width/2+3:width-1])
# # plt.show()
# leftSum = magnitude[0:height/2-3,0:width/2-5].sum()
# rightSum = magnitude[0:height/2-3,width/2+3:width-1].sum()
# print leftSum, rightSum
# print 2*np.pi/x_0, 2*np.pi/y_0
# if leftSum >= rightSum:
#     return A, B, x_0, y_0, phi
# else:
#     return A, B, x_0, y_0, phi