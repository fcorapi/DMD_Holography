#gratingGenerator
#This program creates diffraction gratings that will be projected onto the DMD.
#Frank Corapi (fcorapi@uwaterloo.ca)
#June 14th, 2019

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as img
from PIL import Image

#******************USER INPUTS***************************
xDim = 912  #Length of DMD
yDim = 1140 #Width of DMD
dir = 'C:\Users\Franky\Desktop\UofT Summer 2019\Gratings\CalibrationGratings\PNG\ ' #Directory to save png image file
dirBMP = 'C:\Users\Franky\Desktop\UofT Summer 2019\Gratings\CalibrationGratings\BMP\ ' #Directory to save converted bmp image file
dirBMP1bit = 'C:\Users\Franky\Desktop\UofT Summer 2019\Gratings\CalibrationGratings\BMP_8bit\ '
filename = 'CG' #Filename for the grating
ext1 = '.png' #image extension 1
ext2 = '.bmp' #image extension 2
ext3 = '_8bit.bmp'

period = 2 #period of grating in pixels
angle = 0  #grating angle in radians
#NX = 25 #Number of patches in X-direction
#NY = 30 #Number of patches in Y-direction
r = 57 #radius of calibration patch
#*******************************************************


xvals = np.linspace(0,xDim-1,xDim)
yvals = np.linspace(0,yDim-1,yDim)
xMesh, yMesh = np.meshgrid(xvals,yvals)

def grating(x,y,p,theta):
    g = np.round(0.5*(1+np.cos((2*np.pi/p)*(x*np.cos(theta)+y*np.sin(theta)))))
    return g


for xPos in np.arange(r,xDim+r,2*r):
    for yPos in np.arange(r,yDim+r, 2*r):
        filenameNumber = '_X_'+str(xPos)+'_Y_'+str(yPos)
        g1 = grating(xMesh, yMesh, period, angle)
        g2 = grating(xMesh, yMesh, period, angle)
        for loopx in range(0,xDim-1):
            for loopy in range(0,yDim-1):
                if (loopx - xPos)**2 +(loopy-yPos)**2 > r**2:
                    g1[loopy, loopx] = 0
                if (loopx - ((xDim-2*r)/2))**2 + (loopy-((yDim-2*r)/2))**2 > r**2:
                    g2[loopy,loopx] = 0
        #plt.imshow(g, cmap='gray')
        #plt.show()
        g = g1+g2
        img.imsave(dir+filename+filenameNumber+ext1, g, format = 'png', cmap='gray')

        pngImage = Image.open(dir+filename+filenameNumber+ext1)
        pngImage.save(dirBMP+filename+filenameNumber+ext2)

        oldBmp = Image.open(dirBMP+filename+filenameNumber+ext2)
        newBmp = oldBmp.convert('L')
        newBmp.save(dirBMP1bit+filename+filenameNumber+ext3)
        print 'X = ', xPos, 'Y = ', yPos, ' completed.'
