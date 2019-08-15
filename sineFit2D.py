#sineFit2D
#This program fits a 2-dimensional sine function to an image. Used for determining the phase parameters
#for the DMD interference patterns in order to obtain a calibration phase map. It will also generate an amplitude mad
#and a hologram used for beam shaping.
#Modified the 2D Gaussian Fit code from Scipy Cookbook (https://scipy-cookbook.readthedocs.io/items/FittingData.html)
#Frank Corapi (fcorapi@uwaterloo.ca)
#August 15th, 2019

#Import Directories
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.image as img
from scipy import optimize, interpolate
import csv


#****************************FUNCTIONS*****************************************************************8
#Define the fitting function, a 2-dimensional sine function
def sine2D(A, B, x_0, y_0, phi, wx, wy, gx_0, gy_0):
    A = float(A)
    B = float(B)
    x_0 = float(x_0)
    y_0 = float(y_0)
    phi = float(phi)
    return lambda x,y: (A**2 + B**2 + 2*A*B*np.cos(x_0*x + y_0*y + phi))*np.exp(-np.power(((x-gx_0)/wx),2)-np.power(((y-gy_0)/wy),2))

#Determines the indices of the local minima within a list
def local_min(list):
    minIndices = [i for i,y in enumerate(list)
                  if ((i == 0) or (list[i-1] >= y))
                  and ((i == len(list)-1) or (y < list[i+1]))]
    return minIndices

#This function is used to make an accurate initial guess for the fitting paramters for the 2D sine function
def guessParams(data,xPos, yPos):
    xPos = float(xPos)
    yPos = float(yPos)
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

    #Position of the Center patch (X=1 Y=1 is top left corner of DMD) NEEDS TO CHANGE IF CENTER PATCH CHANGES
    xCen = 8.0
    yCen = 5.0

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
    print deltaX, deltaY, theta,Xi, xPos, yPos
    periodGuessX = (wavelength/cameraPixelSize)/(np.sin(Xi)*np.cos(theta))
    periodGuessY = (wavelength/cameraPixelSize)/(np.sin(Xi)*np.sin(theta))

    x_0 = 2*np.pi/periodGuessX
    y_0 = 2*np.pi/periodGuessY

    phi = 0

    wx = 100
    wy = 100
    gx_0 = 200 #Half the Size of the length of the cropped region
    gy_0 = 175 #Half the Size of the width of the cropped region

    return A, B, x_0, y_0, phi,wx, wy, gx_0, gy_0

#Using the initial guess parameters from guessParams, the optimal parameters are determined using the following function
#A least Squares method is used.
def fitSine2D(data, xPos, yPos):
    params = guessParams(data, xPos, yPos)
    paramBounds = ([0,0,-np.inf,-np.inf,0,0,0,0,0],[np.inf,np.inf,np.inf,np.inf,2*np.pi,np.inf,np.inf,np.inf,np.inf])
    errorfunction = lambda p: np.ravel(sine2D(*p)(*np.indices(data.shape)) -
                                 data)
    p = optimize.least_squares(errorfunction, params, bounds = paramBounds)
    # print success
    return p.x

#Generates a hologram by modulating the phase and the amplitude using the phase difference
#  and amplitude coefficient maps
def hologram(x,y,p,theta,phi,threshold):
    h = 0.5*(1+np.cos((2*np.pi/p)*(x*np.cos(theta)+y*np.sin(theta)) + phi))

    for xLoop in range(np.shape(threshold)[0]):
        for yLoop in range(np.shape(threshold)[1]):
            if h[xLoop,yLoop] < threshold[xLoop, yLoop]:
                h[xLoop,yLoop] = 0
            elif h[xLoop,yLoop] >= threshold[xLoop, yLoop]:
                h[xLoop, yLoop] = 1
    return h

#Just a Gaussian Function
def gaussian(a, x_0, y_0, wx, wy):
    wx = float(wx)
    wy = float(wy)
    return lambda x,y: a*np.exp(-np.power(((x-x_0)/wx),2)-np.power(((y-y_0)/wy),2))

#****************END OF FUNCTIONS**************************************************************************************

#Image Directories
dir = 'C:\Users\Franky\Desktop\UofT Summer 2019\CalibrationImages3 (July 25)\\'
targetDir = 'C:\Users\Franky\Desktop\UofT Summer 2019\Images\\'
targetFilename = 'smallTriforce2'
cropDir = dir + 'Cropped2\\'
fitDir = dir + 'Fitted2\\'
hologramDir = dir + 'Hologram\\'
hologramFilename = 'TestHologram1'
filename = 'CI3_X4Y7'
csvFilename = 'CI3_Params2.csv'
filenamePrefix = 'CI3_'
ext1 = '.bmp'
ext2 = '.png'

#Select the cropping region
cropCoords = (500,300,850,700)

#DMD Image Parameters
xMax = 16 #Maximum X-value in calibration images
yMax = 10 #Maximum Y-value in calibration images
xDim = 912 #Length of DMD image
yDim = 1140#Width of DMD image
xDMDDim = 1290
yDMDDim = 806
#These radii need to match the radii of the calibration patches made using gratingGenerator.py
r = 1140/(2*np.sqrt(2)*10) #radius of calibration patch
r1 = r/np.sqrt(2)
r2 = np.sqrt(2)*r
period = 4 #period of grating in pixels (match original grating)
angle = 0.2 #grating angle in radians (match original grating)

#Program Sequence Initialization
fittingCheck = 0 #set to 1 to fit interference patterns
mapCheck = 1 #set to 1 to generate amplitude and phase maps
hologramCheck = 1 #set to 1 to generate hologram (needs map check set to 1)

#**********************************FITTING CALIBRATION IMAGES*****************************************************
if fittingCheck == 1:
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
            print "Params:", params
            # print params
            # print 2 * np.pi / params[2], 2 * np.pi / params[3]
            fit = sine2D(*params)
            plt.figure()
            plt.imshow(croppedData, cmap = 'viridis')
            plt.colorbar()
            plt.contour(fit(*np.indices(croppedData.shape)), cmap=plt.cm.get_cmap('copper'), linewidths = 2)
            plt.savefig(fitDir + filenamePrefix + 'X' + str(xLoop) + 'Y' + str(yLoop) + '_Fit' + ext2)
            plt.close()
            headers = ['A', 'B', 'x_0', 'y_0', 'phi']
            if xLoop == 1 and yLoop == 1:
                # with open(dir+csvFilename, "wb") as output:
                #     writer = csv.writer(output, lineterminator = '\n')
                #     writer.writerow(headers)
                with open(dir + csvFilename, "wb") as output:
                    writer = csv.writer(output, lineterminator='\n')
                    writer.writerow(params)
            else:
                with open(dir+csvFilename, "ab") as output:
                    writer = csv.writer(output, lineterminator = '\n')
                    writer.writerow(params)

#*************************************************GENRATING PHASE AND AMPLITUDE MAPS*******************************
if mapCheck == 1:
    AList = []
    BList = []
    phiList = []
    with open(dir+csvFilename, 'r') as csvFile:
        reader = csv.reader(csvFile)
        for row in reader:
            AList.append(float(row[0]))
            BList.append(float(row[1]))
            phiList.append(float(row[4]))

    xValues = np.arange(r1,xDim+r1-1,2*r1)
    yValues = np.arange(r2, yDim+r2-1,2*r2)
    xMesh, yMesh = np.meshgrid(xValues, yValues)

    phiMap = np.zeros([len(yValues), len(xValues)])
    ampMap = np.zeros([len(yValues), len(xValues)])

    listLoop = 0
    for xLoop in range(0,xMax):
        for yLoop in range(0,yMax):
            if xLoop == 7 and yLoop == 4:
                phiMap[yLoop,xLoop] = 0
            else:
                phiMap[yLoop,xLoop] = phiList[listLoop]
                listLoop = listLoop + 1

    listLoop = 0
    meanAList = []
    for xLoop in range(0, xMax):
        for yLoop in range(0, yMax):
            if xLoop == 7 and yLoop == 4:
                ampMap[yLoop, xLoop] = 0
            else:
                ampMap[yLoop, xLoop] = AList[listLoop]*BList[listLoop]
                meanAList.append(max(AList[listLoop],BList[listLoop]))
                listLoop = listLoop + 1

    #Vertical Phase Unwrapping
    unwrapPhiMap = np.zeros(np.shape(phiMap))
    unwrapPhiMap[:,:] = phiMap
    for wrapLoopx in range(0,xMax):
        for wrapLoopy in range(2,yMax):
            slope = unwrapPhiMap[wrapLoopy-1,wrapLoopx] - unwrapPhiMap[wrapLoopy-2,wrapLoopx]
            expectedPhi = unwrapPhiMap[wrapLoopy-1,wrapLoopx] + slope
            measuredPhi = unwrapPhiMap[wrapLoopy,wrapLoopx]
            m = round((expectedPhi-measuredPhi)/(2*np.pi))
            unwrapPhiMap[wrapLoopy, wrapLoopx] = unwrapPhiMap[wrapLoopy,wrapLoopx]+2*np.pi*m

    #Horizontal Phase Unwrapping
    unwrapPhiMap2 = np.zeros(np.shape(unwrapPhiMap))
    unwrapPhiMap2[:, :] = unwrapPhiMap
    for wrapLoopx in range(2, xMax):
        for wrapLoopy in range(0, yMax):
            slope = unwrapPhiMap2[wrapLoopy, wrapLoopx-1] - unwrapPhiMap2[wrapLoopy, wrapLoopx-2]
            expectedPhi = unwrapPhiMap2[wrapLoopy, wrapLoopx-1] + slope
            measuredPhi = unwrapPhiMap2[wrapLoopy, wrapLoopx]
            m = round((expectedPhi - measuredPhi) / (2 * np.pi))
            unwrapPhiMap2[wrapLoopy, wrapLoopx] = unwrapPhiMap2[wrapLoopy, wrapLoopx] + 2 * np.pi * m

    #Interpolation of Phase Map
    interpPhiMapObject = interpolate.interp2d(xValues, yValues, unwrapPhiMap2, kind='cubic')

    newXValuesPhi = np.arange(r1, xDim - r1, 1)
    newYValuesPhi = np.arange(r2, yDim - r2, 1)

    interpPhiMap = interpPhiMapObject(newXValuesPhi, newYValuesPhi)

    interpPhiMapResized = np.zeros([yDim, xDim])
    interpPhiMapResized[int(r2):interpPhiMap.shape[0] + int(r2), int(r1):interpPhiMap.shape[1] + int(r1)] = interpPhiMap

    #Determines an approximation for the amplitude of the center calibration patch
    meanA = np.average(meanAList)
    ampMap[4,7] = meanA**2

    #Interpolation of Amplitude Map
    interpAmpMapObject = interpolate.interp2d(xValues,yValues,ampMap, kind='cubic')

    newXValues = np.arange(r1, xDim-r1, 1)
    newYValues = np.arange(r2, yDim-r2, 1)

    interpAmpMap = interpAmpMapObject(newXValues, newYValues)

    normAmpMap = np.zeros(np.shape(interpAmpMap))
    normAmpMap[:,:] = interpAmpMap/np.amax(interpAmpMap)

    interpAmpMapResized = np.zeros([yDim, xDim])
    interpAmpMapResized[int(r2):interpAmpMap.shape[0]+int(r2),int(r1):interpAmpMap.shape[1]+int(r1)] = interpAmpMap

    normAmpMapResized = np.zeros(np.shape(interpAmpMapResized))
    normAmpMapResized[:,:] = interpAmpMapResized/np.amax(interpAmpMapResized)

    # #PLOTTING
    # plt.figure(1)
    # plt.imshow(phiMap,extent=(0, xDim, yDim, 0), interpolation='none', cmap='rainbow')
    # plt.title('Phase')
    # plt.xlabel('X')
    # plt.ylabel('Y')
    # plt.colorbar()
    #
    # plt.figure(2)
    # plt.imshow(unwrapPhiMap, extent=(0, xDim, yDim, 0), interpolation='none', cmap='rainbow')
    # plt.title('Phase, Unwrapped Vertically')
    # plt.xlabel('X')
    # plt.ylabel('Y')
    # plt.colorbar()
    #
    # plt.figure(3)
    # plt.imshow(unwrapPhiMap2, extent=(0, xDim, yDim, 0), interpolation='none', cmap='rainbow')
    # plt.title('Phase, Unwrapped')
    # plt.xlabel('X')
    # plt.ylabel('Y')
    # plt.colorbar()
    #
    # plt.figure(4)
    # plt.imshow(interpPhiMapResized, extent=(0, xDim, yDim, 0), interpolation='none', cmap='rainbow')
    # plt.title('Phase, Unwrapped (Interpolated)')
    # plt.xlabel('X')
    # plt.ylabel('Y')
    # plt.colorbar()
    #
    # plt.figure(5)
    # plt.imshow(interpAmpMapResized, extent=(0, xDim, yDim, 0), interpolation='none',
    #            cmap='rainbow')
    # #origin = 'upper', if using contourf
    # plt.title('Amplitude')
    # plt.xlabel('X')
    # plt.ylabel('Y')
    # #plt.gca().invert_yaxis() if using contourf
    # plt.colorbar()
    # plt.show()

#******************************HOLOGRAM GENERATION*****************************************************
if hologramCheck == 1:

    #If wanting to create a gaussian image, use this code
    # xVals = np.arange(0, xDim, 1)
    # yVals = np.arange(0, yDim, 1)
    # xMesh, yMesh = np.meshgrid(xVals, yVals)
    #
    # gauss = gaussian(1, xDim/2, yDim/2, 1, 1)(xMesh, yMesh)

    #Obtain Target Image
    targetImage = Image.open(targetDir+targetFilename+ext1)
    target = np.array(targetImage) #Use this or gauss for testing

    #Calculate the amplitude and phase of the FT of the Target Image
    targetFT = np.fft.fft2(target)
    targetFTShift = np.fft.fftshift(targetFT)
    targetAmp = np.abs(targetFTShift)/np.amax(np.abs(targetFTShift))
    targetPhase = np.arctan(np.imag(targetFTShift)/(np.real(targetFTShift)+1e-10))

    print np.shape(targetFTShift)

    #Determine the amplitude coefficient map needed for amplitude modulation
    ampCoeff = targetAmp/normAmpMap
    ampCoeffNorm = ampCoeff/np.amax(ampCoeff)
    # print np.shape(ampCoeff)

    #Resize the map by padding outside with zeros
    ampCoeffNormResized = np.zeros([yDim, xDim])
    ampCoeffNormResized[int(r2):ampCoeffNorm.shape[0] + int(r2), int(r1):ampCoeffNorm.shape[1] + int(r1)] = ampCoeffNorm

    phaseDifference = targetPhase - interpPhiMap

    phaseDifferenceResized = np.zeros([yDim, xDim])
    phaseDifferenceResized[int(r2):phaseDifference.shape[0]+int(r2),int(r1):phaseDifference.shape[1]+int(r1)] = phaseDifference

    ampThreshold = np.zeros(np.shape(ampCoeffNorm))
    ampThreshold[:,:] = 1-ampCoeffNorm/2

    hologramXVals = np.arange(r1, xDim-r1, 1)
    hologramYVals = np.arange(r2, yDim-r2, 1)
    hologramXMesh, hologramYMesh = np.meshgrid(hologramXVals, hologramYVals)

    h = hologram(hologramXMesh, hologramYMesh, period, angle, phaseDifference, ampThreshold)

    holo = np.zeros([yDim,xDim])
    holo[int(r2):h.shape[0]+int(r2),int(r1):h.shape[1]+int(r1)] = h

    hologramImage = Image.new('1', (xDim,yDim))
    pixelImage = hologramImage.load()
    for i in range(hologramImage.size[0]):
        for j in range(hologramImage.size[1]):
            pixelImage[i, j] = (holo[j][i],)

    hologramImage.save(hologramDir+hologramFilename+ext1)

    #Plot Target Image
    plt.figure(1)
    plt.imshow(target)
    plt.colorbar()

    # Plot Amplitude and Phase of Target Image
    plt.figure(2)
    plt.imshow(targetAmp)
    plt.colorbar()

    plt.figure(3)
    plt.imshow(targetPhase)
    plt.colorbar()

    plt.figure(4)
    plt.imshow(phaseDifferenceResized)
    plt.colorbar()

    plt.figure(5)
    plt.imshow(ampCoeffNormResized)
    plt.colorbar()
    # plt.clim(0, 0.02)
    plt.show()



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