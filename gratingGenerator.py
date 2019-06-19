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
dir = 'C:\Users\Franky\Desktop\UofT Summer 2019\Gratings\ ' #Directory to save png image file
dirBMP = 'C:\Users\Franky\Desktop\UofT Summer 2019\Gratings\BMP\ ' #Directory to save converted bmp image file
filename = 'grating1' #Filename for the grating
ext1 = '.png' #image extension 1
ext2 = '.bmp' #image extension 2
ext3 = '_1bit.bmp'

period = 10 #period of grating in pixels
angle = 0  #grating angle in radians
#*******************************************************


xvals = np.linspace(0,xDim-1,xDim)
yvals = np.linspace(0,yDim-1,yDim)
xMesh, yMesh = np.meshgrid(xvals,yvals)

def grating(x,y,p,theta):
    g = np.round(0.5*(1+np.cos((2*np.pi/p)*(x*np.cos(theta)+y*np.sin(theta)))))
    return g

gTest = grating(xMesh, yMesh,period,angle)
print np.shape(xMesh)
print gTest
plt.imshow(gTest, cmap='gray')
plt.show()

img.imsave(dir+filename+ext1, gTest, format = 'png', cmap='gray')

pngImage = Image.open(dir+filename+ext1)
pngImage.save(dirBMP+filename+ext2)

oldBmp = Image.open(dirBMP+filename+ext2)
newBmp = oldBmp.convert('1')
newBmp.save(dirBMP+filename+ext3)
