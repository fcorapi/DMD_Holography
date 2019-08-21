#fftImage
#This program takes the fourier transform of an image. It is simply used for testing before implementation in main
#program.
#Frank Corapi (fcorapi@uwaterloo.ca)
#Last Modified: August 20th, 2019

import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.image as img

def gaussian(a, x_0, y_0, wx, wy):
    wx = float(wx)
    wy = float(wy)
    return lambda x,y: a*np.exp(-np.power(((x-x_0)/wx),2)-np.power(((y-y_0)/wy),2))

xVals = np.arange(-1,1,0.001)
yVals = np.arange(-1,1,0.001)
xMesh,yMesh = np.meshgrid(xVals,yVals)

gauss = gaussian(1,0,0,0.01,0.01)(xMesh,yMesh)
# plt.figure(1)
# plt.imshow(gauss, extent= (xVals.min(),xVals.max(),yVals.min(),yVals.max()))
# plt.colorbar()



dir = 'C:\Users\Franky\Desktop\UofT Summer 2019\Images\\'
filename = 'startrek'
filename2 = 'Gaussian2_Cropped_FT'
ext1 = '.bmp'
ext2 = '.bmp'
ext3 = '.png'
ext4 = '.jpg'
xDim = 912  #Length of DMD
yDim = 1140 #Width of DMD

initialImg = Image.open(dir+filename+ext1)
np_im = np.array(initialImg)
print np_im
# npImage = np.asarray(initialImg)
# print npImage
# initialImg.save(dir+filename+ext2)
# bmpImage = Image.open(dir+filename+ext2)
# bmpImg = bmpImage.convert('L')
# bmpImg.save(dir+filename+ext2)
plt.figure(1)
plt.imshow(np_im)
plt.colorbar()


ft = np.fft.fft2(np_im)
ftShift = np.fft.fftshift(ft)
magnitude = np.abs(ftShift)/np.amax(np.abs(ftShift))
magnitude2 = np.fft.fftshift(magnitude)
phase = np.arctan(np.imag(ftShift)/np.real(ftShift))
phase2 = np.fft.fftshift(phase)

# img.imsave(dir+filename2+ext3, magnitude2,  format = 'png')
# ftpng = Image.open(dir+filename2+ext3)
# ftpng.save(dir+filename2+ext2)
# bmpFT = Image.open(dir+filename2+ext2)
# bmpFT = bmpFT.convert('L')
# bmpFT.save(dir+filename2+ext2)

plt.figure(2)
plt.contourf(magnitude,100)
plt.colorbar()

plt.figure(3)
plt.contourf(phase,10)
plt.colorbar()

# ift = np.fft.ifftshift(ftShift)
# ift = np.fft.ifft2(ift)
# plt.figure(4)
# plt.imshow(np.abs(ift))
# plt.colorbar()
plt.show()
