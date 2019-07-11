#gratingGenerator
#This program creates diffraction gratings that will be projected onto the DMD.
#Frank Corapi (fcorapi@uwaterloo.ca)
#June 20th, 2019

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as img
from PIL import Image

#******************USER INPUTS***************************
xDim = 912 #Length of DMD
yDim = 1140#Width of DMD
dir = 'C:\Users\Franky\Desktop\UofT Summer 2019\Gratings\CalibrationGratings\PNG\\' #Directory to save png image file
dirBMP = 'C:\Users\Franky\Desktop\UofT Summer 2019\Gratings\CalibrationGratings\BMP\\' #Directory to save converted bmp image file
dirBMP8bit = 'C:\Users\Franky\Desktop\UofT Summer 2019\Gratings\CalibrationGratings\BMP_8bit\\'
dirBMP1bit = 'C:\Users\Franky\Desktop\UofT Summer 2019\Gratings\CalibrationGratings4\\'
filename = 'CG4' #Filename for the grating
ext1 = '.png' #image extension 1
ext2 = '.bmp' #image extension 2
ext3 = '_8bit.bmp' #image extension 2
ext4 = '_1bit.bmp' #image extension 2

period = 4 #period of grating in pixels
angle = 0.2 #grating angle in radians
#NX = 25 #Number of patches in X-direction
#NY = 30 #Number of patches in Y-direction
r = 1140/(2*np.sqrt(2)*5) #radius of calibration patch
r1 = r/np.sqrt(2)
r2 = np.sqrt(2)*r
#*******************************************************


xvals = np.linspace(0,xDim-1,xDim)
yvals = np.linspace(0,yDim-1,yDim)
xMesh, yMesh = np.meshgrid(xvals,yvals)

def grating(x,y,p,theta):
    g = np.round(0.5*(1+np.cos((2*np.pi/p)*(x*np.cos(theta)+y*np.sin(theta)))))
    return g

def truncate(n, decimals):
    multiplier = 10.0 ** decimals
    return float(int(n * multiplier) / multiplier)

for xPos in np.arange(r1,xDim+r1-1,2*r1):
    for yPos in np.arange(r2,yDim+r2-1, 2*r2):
        filenameNumber = '_X_'+str(truncate(xPos,2))+'_Y_'+str(truncate(yPos,2))
        g1 = grating(xMesh, yMesh, period, angle)
        g2 = grating(xMesh, yMesh, period, angle)
        for loopx in range(0,xDim):
            for loopy in range(0,yDim):
                if ((loopx - xPos)/r1)**2 +((loopy-yPos)/r2)**2 > 1:
                    g1[loopy, loopx] = 0
                if ((loopx - ((xDim-2*r1)/2))/r1)**2 + ((loopy-((((yDim-2*r2)/2))+r2))/r2)**2 > 1:
                    g2[loopy,loopx] = 0
        g = g1+g2
        directImage = Image.new('1', (xDim,yDim))
        pixelImage = directImage.load()
        for i in range(directImage.size[0]):
            for j in range(directImage.size[1]):
                pixelImage[i,j] = (g[j][i],)
        directImage.save(dirBMP1bit+filename+filenameNumber+ext4)

        print 'X = ', truncate(xPos,2), 'Y = ', truncate(yPos,2), ' completed.'

# g = grating(xMesh,yMesh,period,angle)
# directImage = Image.new('1', (xDim,yDim))
# pixelImage = directImage.load()
# for i in range(directImage.size[0]):
#     for j in range(directImage.size[1]):
#         pixelImage[i, j] = (g[j][i],)
#         # if (j > 200 and j<624) and (i >200 and i <412):
#         #     pixelImage[i,j] = (1,)
#         # else:
#         #     pixelImage[i, j] = (0,)
directImage.save(dirBMP1bit+filename+ext4)
#************************OLD ATTEMPTS**********************************
        #plt.imshow(g, cmap='gray')
        #plt.show()

        #img.imsave(dir+filename+filenameNumber+ext1, g, format = 'png', cmap='gray')
        # directImage = Image.fromarray(g,'L')
        # directImage.show()

        #pngImage = Image.open(dir+filename+filenameNumber+ext1)
        #pngImage.save(dirBMP+filename+filenameNumber+ext2)

        #oldBmp = Image.open(dirBMP+filename+filenameNumber+ext2)
        #newBmp = oldBmp.convert('L')
        #newBmp.save(dirBMP8bit+filename+filenameNumber+ext3)
        #newerBmp = Image.open(dirBMP8bit+filename+filenameNumber+ext3)
        #newestBmp = newerBmp.convert('1')
        #newestBmp.save(dirBMP1bit+filename+filenameNumber+ext4)