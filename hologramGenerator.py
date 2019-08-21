#hologramGenerator
#This program fits a 2-dimensional sine function, multiplied by a 2D gaussian to an image. Used for determining the phase parameters
#for the DMD interference patterns in order to obtain a calibration phase map. It then generates an amplitude map
#and a hologram used for beam shaping.
#Modified the 2D Gaussian Fit code from Scipy Cookbook (https://scipy-cookbook.readthedocs.io/items/FittingData.html)
#Frank Corapi (fcorapi@uwaterloo.ca)
#Last Updated: August 21st, 2019

#Import Directories
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.image as img
from scipy import optimize, interpolate
import csv


#****************************FUNCTIONS*****************************************************************8
#Define the fitting function, a 2-dimensional sine function multipled by a 2-dimensional gaussian function
def sine2D(A, B, x_0, y_0, phi, wx, wy, gx_0, gy_0):
    A = float(A)
    B = float(B)
    x_0 = float(x_0)
    y_0 = float(y_0)
    phi = float(phi)
    return lambda x,y: (A**2 + B**2 + 2*A*B*np.cos(x_0*x + y_0*y + phi))*np.exp(-np.power(((x-gx_0)/wx),2)-np.power(((y-gy_0)/wy),2))

#Determines the indices of the local minima within a list (Not used anymore)
def local_min(list):
    minIndices = [i for i,y in enumerate(list)
                  if ((i == 0) or (list[i-1] >= y))
                  and ((i == len(list)-1) or (y < list[i+1]))]
    return minIndices

#This function is used to make an accurate initial guess for the fitting paramters for the 2D sine function.
#It uses the relative positions of the calibration patches to calculate a guess for the period and angle of the
#interference pattern. The guess for the amplitude assumes that A and B are equal and that the max value of the image
#equal to (A+B)^2.
def guessParams(data,xPos, yPos, xCen, yCen):
    xPos = float(xPos)
    yPos = float(yPos)
    xCen = float(xCen)
    yCen = float(yCen)
    amplitude = data.max()

    A = np.sqrt(amplitude)/2
    B = A

    #Optical System Parameters (NEEDS TO BE MODIFIED IF PATCHES OR SYSTEM CHANGES)
    wavelength = 0.760 #microns
    dPatch = 0.0616 #cm (distance between two adjacent patch centers)
    cameraPixelSize = 5.2 #microns/pixel
    f1 = 10.0 #cm (focal length of first lens after DMD)
    f2 = 4.0 #cm (focal length of second lens after DMD)
    f3 = 30.0 #cm (focal length of third lens after DMD)

    #Position of the Center patch (X=1 Y=1 is top left corner of DMD) NEEDS TO CHANGE IF CENTER PATCH CHANGES
    #MUST BE CONSISTENT WITH gratingGenerator CENTER PATCH LOCATION
    #NO LONGER UPDATE HERE BUT UPDATED AFTER FUNCTION DEFINITIONS (LINES 161 AND 162)
    # xCen = 8.0 #must be float
    # yCen = 5.0 #must be float

    deltaX = xCen - xPos
    if deltaX == 0:
        deltaX = 1e-3

    deltaY = yCen - yPos
    if deltaY == 0:
        deltaY = 1e-3

    #Determines which quandrant the calibration patch resides in
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
    # print deltaX, deltaY, theta, Xi, xPos, yPos
    periodGuessX = (wavelength/cameraPixelSize)/(np.sin(Xi)*np.cos(theta))
    periodGuessY = (wavelength/cameraPixelSize)/(np.sin(Xi)*np.sin(theta))

    x_0 = 2*np.pi/periodGuessX
    y_0 = 2*np.pi/periodGuessY

    phi = 0

    wx = 100 #Tweaked by just guessing until gaussian fit became accurate
    wy = 100 #Tweaked by just guessing until gaussian fit became accurate
    gx_0 = 200 #Half the Size of the length of the cropped region
    gy_0 = 175 #Half the Size of the width of the cropped region

    # Returns the guess for the parameters
    return A, B, x_0, y_0, phi,wx, wy, gx_0, gy_0

#Using the initial guess parameters from guessParams, the optimal parameters are determined using the following function
#A least Squares method is used.
def fitSine2D(data, xPos, yPos, xCen, yCen):
    #Obtain guess parameters
    params = guessParams(data, xPos, yPos, xCen, yCen)
    #Set the bounds on what the determined fit parameters can be
    paramBounds = ([0,0,-np.inf,-np.inf,0,0,0,0,0],[np.inf,np.inf,np.inf,np.inf,2*np.pi,np.inf,np.inf,np.inf,np.inf])
    #Define the error function to be minimized (difference between fit function and data)
    errorfunction = lambda p: np.ravel(sine2D(*p)(*np.indices(data.shape)) -
                                 data)
    #Minimize error function
    p = optimize.least_squares(errorfunction, params, bounds = paramBounds)
    #Returns the fit parameters
    return p.x

#Generates a hologram by modulating the phase and the amplitude using the phase difference
#and amplitude coefficient maps.
def hologram(x,y,p,theta,phi,threshold):
    h = 0.5*(1+np.cos((2*np.pi/p)*(x*np.cos(theta)+y*np.sin(theta)) + phi)) #Grating Equation

    #Check whether the value of the grating equation at each pixel is greater or less than the amplitude threshold
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

#*********************************END OF FUNCTIONS************************************************************

#Image Directories
dir = 'C:\Users\Franky\Desktop\UofT Summer 2019\CalibrationImages3 (July 25)\\' #Directory containing calibration imgs
targetDir = 'C:\Users\Franky\Desktop\UofT Summer 2019\Images\\' #Directory containing target image
#*********
#NOTE: The target image should be the size of the interpolated phase and amplitude maps, NOT THE DMD IMAGE SIZE.
#*********
targetFilename = 'triforce_RS2' #Filename of target excluding extensions (USER INPUT REQUIRED)
cropDir = dir + 'Cropped2\\' #Directory where cropped images will be saved
fitDir = dir + 'Fitted2\\' #Driectory where fitted images will be saved
hologramDir = dir + 'Hologram\\' #Directory where hologram grating will be saved
hologramFilename = 'TestHologram2' #Filename of hologram grating exluding extension (USER INPUT REQUIRED)
csvFilename = 'CI3_Params2.csv' #Set the filename of the CSV file containing the fit parameters (USER INPUT REQUIRED)
filenamePrefix = 'CI3_' #Used to help read in the calibration images
ext1 = '.bmp' #filename extension
ext2 = '.png' #filename extension

#Select the cropping region
cropCoords = (500,300,850,700)

#DMD Image Parameters
xMax = 16 #Maximum X-value in calibration images (Number of patches in X direction) (keep as int)
yMax = 10 #Maximum Y-value in calibration images (Number of patches in Y direction) (keep as int)
xCenter = 8 #X-Location of Center Patch (keep as int)
yCenter = 5 #Y-Location of Center Patch (keep as int)
xDim = 912 #Length of DMD image (keep as int)
yDim = 1140#Width of DMD image (keep as int)
xDMDDim = 1290 #Not used
yDMDDim = 806 #Not used
#These radii need to match the radii of the calibration patches made using gratingGenerator.py
#****************MUST BE CONSISTENT WITH gratingGenerator**************************
r = 1140/(2*np.sqrt(2)*10) #radius of calibration patch
r1 = r/np.sqrt(2) #Scaling of radius due to DMD mirror orientations
r2 = np.sqrt(2)*r #Scaling of radius due to DMD mirror orientations
period = 4 #period of grating in pixels (match original grating)
angle = 0.2 #grating angle in radians (match original grating)
#*********************************************************************************

#Program Sequence Initialization
fittingCheck = 0 #set to 1 to fit interference patterns
mapCheck = 1 #set to 1 to generate amplitude and phase maps
hologramCheck = 1 #set to 1 to generate hologram (needs map check set to 1)

#**********************************FITTING CALIBRATION IMAGES*****************************************************
if fittingCheck == 1:
    for xLoop in range(1,xMax+1):
        for yLoop in range(1,yMax+1):
            if xLoop == xCenter and yLoop == yCenter: #Skips the center patch (May be modified in the future)
                continue
            image = Image.open(dir + filenamePrefix + 'X' + str(xLoop) + 'Y' + str(yLoop) + ext1)
            # data = np.array(image)[:, :, 0]
            # plt.imshow(data)
            # plt.show()

            #Crop and Save Image
            croppedImage = image.crop(cropCoords)
            croppedImage.save(cropDir + filenamePrefix + 'X' + str(xLoop) + 'Y' + str(yLoop) + '_Cropped'+ ext1)
            croppedData = np.array(croppedImage)[:, :, 0]

            params = fitSine2D(croppedData, xLoop, yLoop, xCenter, yCenter)
            print "Params:", params
            # print params
            # print 2 * np.pi / params[2], 2 * np.pi / params[3]
            fit = sine2D(*params)

            #***Save Fit Image***
            plt.figure()
            plt.imshow(croppedData, cmap = 'viridis')
            plt.colorbar()
            plt.contour(fit(*np.indices(croppedData.shape)), cmap=plt.cm.get_cmap('copper'), linewidths = 2)
            plt.savefig(fitDir + filenamePrefix + 'X' + str(xLoop) + 'Y' + str(yLoop) + '_Fit' + ext2)
            plt.close()
            #*******************

            # headers = ['A', 'B', 'x_0', 'y_0', 'phi'] #don't use anymore
            #Writes parameters for each fit into a csv file
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
    #Read CSV file in
    with open(dir+csvFilename, 'r') as csvFile:
        reader = csv.reader(csvFile)
        for row in reader:
            AList.append(float(row[0]))
            BList.append(float(row[1]))
            phiList.append(float(row[4]))

    xValues = np.arange(r1,xDim+r1-1,2*r1)
    yValues = np.arange(r2, yDim+r2-1,2*r2)
    xMesh, yMesh = np.meshgrid(xValues, yValues)

    #Create empty phase and amplitude map
    phiMap = np.zeros([len(yValues), len(xValues)])
    ampMap = np.zeros([len(yValues), len(xValues)])

    #Fills phi map with phi parameter from fit
    listLoop = 0
    for xLoop in range(0,xMax):
        for yLoop in range(0,yMax):
            if xLoop == xCenter-1 and yLoop == yCenter-1:
                phiMap[yLoop,xLoop] = 0
            else:
                phiMap[yLoop,xLoop] = phiList[listLoop]
                listLoop = listLoop + 1

    #Fills amplitude map with the product of the A and B parameter from fit
    listLoop = 0
    meanAList = []
    for xLoop in range(0, xMax):
        for yLoop in range(0, yMax):
            if xLoop == xCenter-1 and yLoop == yCenter-1:
                ampMap[yLoop, xLoop] = 0
            else:
                ampMap[yLoop, xLoop] = AList[listLoop]*BList[listLoop]
                meanAList.append(max(AList[listLoop],BList[listLoop])) #Keeps track of the intensity value from the center patch
                listLoop = listLoop + 1

    #*********************Vertical Phase Unwrapping*********************
    unwrapPhiMap = np.zeros(np.shape(phiMap))
    unwrapPhiMap[:,:] = phiMap
    for wrapLoopx in range(0,xMax):
        for wrapLoopy in range(2,yMax):
            slope = unwrapPhiMap[wrapLoopy-1,wrapLoopx] - unwrapPhiMap[wrapLoopy-2,wrapLoopx]
            expectedPhi = unwrapPhiMap[wrapLoopy-1,wrapLoopx] + slope
            measuredPhi = unwrapPhiMap[wrapLoopy,wrapLoopx]
            m = round((expectedPhi-measuredPhi)/(2*np.pi))
            unwrapPhiMap[wrapLoopy, wrapLoopx] = unwrapPhiMap[wrapLoopy,wrapLoopx]+2*np.pi*m
    #*******************************************************************

    #*************Horizontal Phase Unwrapping***************************
    unwrapPhiMap2 = np.zeros(np.shape(unwrapPhiMap))
    unwrapPhiMap2[:, :] = unwrapPhiMap
    for wrapLoopx in range(2, xMax):
        for wrapLoopy in range(0, yMax):
            slope = unwrapPhiMap2[wrapLoopy, wrapLoopx-1] - unwrapPhiMap2[wrapLoopy, wrapLoopx-2]
            expectedPhi = unwrapPhiMap2[wrapLoopy, wrapLoopx-1] + slope
            measuredPhi = unwrapPhiMap2[wrapLoopy, wrapLoopx]
            m = round((expectedPhi - measuredPhi) / (2 * np.pi))
            unwrapPhiMap2[wrapLoopy, wrapLoopx] = unwrapPhiMap2[wrapLoopy, wrapLoopx] + 2 * np.pi * m
    #*******************************************************************

    #***********************Interpolation of Phase Map********************************
    interpPhiMapObject = interpolate.interp2d(xValues, yValues, unwrapPhiMap2, kind='cubic')

    newXValuesPhi = np.arange(r1, xDim - r1, 1)
    newYValuesPhi = np.arange(r2, yDim - r2, 1)

    interpPhiMap = interpPhiMapObject(newXValuesPhi, newYValuesPhi)

    #Resizes the interpolated phi map to be the size of the DMD acceptable image
    interpPhiMapResized = np.zeros([yDim, xDim])
    interpPhiMapResized[int(r2):interpPhiMap.shape[0] + int(r2), int(r1):interpPhiMap.shape[1] + int(r1)] = interpPhiMap
    # **********************************************************************************

    #Determines an approximation for the amplitude of the center calibration patch
    meanA = np.average(meanAList)
    ampMap[xCenter-1, yCenter-1] = meanA**2

    #*****************************Interpolation of Amplitude Map*******************************
    interpAmpMapObject = interpolate.interp2d(xValues,yValues,ampMap, kind='cubic')

    newXValues = np.arange(r1, xDim-r1, 1)
    newYValues = np.arange(r2, yDim-r2, 1)

    interpAmpMap = interpAmpMapObject(newXValues, newYValues)

    #Normalizes the amplitude map
    normAmpMap = np.zeros(np.shape(interpAmpMap))
    normAmpMap[:,:] = interpAmpMap/np.amax(interpAmpMap)

    # Resizes the interpolated amplitude map to be the size of the DMD acceptable image
    interpAmpMapResized = np.zeros([yDim, xDim])
    interpAmpMapResized[int(r2):interpAmpMap.shape[0]+int(r2),int(r1):interpAmpMap.shape[1]+int(r1)] = interpAmpMap

    # Resizes the interpolated and normalized amplitude map to be the size of the DMD acceptable image
    normAmpMapResized = np.zeros(np.shape(interpAmpMapResized))
    normAmpMapResized[:,:] = interpAmpMapResized/np.amax(interpAmpMapResized)
    #*****************************************************************************************************

    # #*********PLOTTING***********
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

    #******************Obtain Target Image***************************************
    targetImage = Image.open(targetDir+targetFilename+ext1)
    target = np.array(targetImage) #Use targetImage or gauss for testing

    #Calculate the amplitude and phase of the FT of the Target Image
    targetFT = np.fft.fft2(target)
    targetFTShift = np.fft.fftshift(targetFT)
    targetAmp = np.abs(targetFTShift)/np.amax(np.abs(targetFTShift)) #NORMALIZED
    targetPhase = np.arctan(np.imag(targetFTShift)/(np.real(targetFTShift)+1e-10))

    print np.shape(targetFTShift) #This should be the size of the interpolated phase and amplitude maps

    #Determine the amplitude coefficient map needed for amplitude modulation
    ampCoeff = targetAmp/normAmpMap #NORMALIZEDTARGET/NORMALIZEDAMPMAP
    ampCoeffNorm = ampCoeff/np.amax(ampCoeff)
    # print np.shape(ampCoeff)

    #*********************Resize the map by padding outside with zeros*********************************
    ampCoeffNormResized = np.zeros([yDim, xDim])
    ampCoeffNormResized[int(r2):ampCoeffNorm.shape[0] + int(r2), int(r1):ampCoeffNorm.shape[1] + int(r1)] = ampCoeffNorm

    #******Determine the phase difference map needed for phase modulation*********
    phaseDifference = targetPhase - interpPhiMap

    # *********************Resize the map by padding outside with zeros*********************************
    phaseDifferenceResized = np.zeros([yDim, xDim])
    phaseDifferenceResized[int(r2):phaseDifference.shape[0]+int(r2),int(r1):phaseDifference.shape[1]+int(r1)] = phaseDifference

    #Calculate the amplitude threshold which will be used to determine which regions of the hologram will have values of zero
    ampThreshold = np.zeros(np.shape(ampCoeffNorm))
    ampThreshold[:,:] = 1-ampCoeffNorm/2

    #********************Generated the hologram************************
    hologramXVals = np.arange(r1, xDim-r1, 1)
    hologramYVals = np.arange(r2, yDim-r2, 1)
    hologramXMesh, hologramYMesh = np.meshgrid(hologramXVals, hologramYVals)
    #Use non-rezied phase difference and ampThreshold
    h = hologram(hologramXMesh, hologramYMesh, period, angle, phaseDifference, ampThreshold)
    #*******************************************************************

    # *********************Resize the hologram by padding outside with zeros*********************************
    holo = np.zeros([yDim,xDim])
    holo[int(r2):h.shape[0]+int(r2),int(r1):h.shape[1]+int(r1)] = h

    #Save Hologram grating to an image
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

    #Plot Target Phase
    plt.figure(3)
    plt.imshow(targetPhase)
    plt.colorbar()

    #Plot Phase Difference Resized
    plt.figure(4)
    plt.imshow(phaseDifferenceResized)
    plt.colorbar()

    #Plot Amplitude Coefficient Map Resized and Normalized
    plt.figure(5)
    plt.imshow(ampCoeffNormResized)
    plt.colorbar()
    # plt.clim(0, 0.02)
    plt.show()
