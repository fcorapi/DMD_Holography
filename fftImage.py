#fftImage
#This program takes the fourier transform of an image.
#Frank Corapi (fcorapi@uwaterloo.ca)
#June 21th, 2019

import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.image as img

dir = 'C:\Users\Franky\Desktop\UofT Summer 2019\ImagePlane\\'
filename = 'Gaussian2_Cropped'
filename2 = 'Gaussian2_Cropped_FT'
ext1 = '.bmp'
ext2 = '.bmp'
ext3 = '.png'
xDim = 912  #Length of DMD
yDim = 1140 #Width of DMD

initialImg = Image.open(dir+filename+ext1)
# initialImg.save(dir+filename+ext2)
# bmpImage = Image.open(dir+filename+ext2)
# bmpImg = bmpImage.convert('L')
# bmpImg.save(dir+filename+ext2)
plt.imshow(initialImg)
plt.show()

ft = np.fft.fft2(initialImg)
ftShift = np.fft.fftshift(ft)
magnitude = 200*np.log(np.abs(ft))
phase = np.arctan(np.imag(ft)/np.real(ft))

img.imsave(dir+filename2+ext3, magnitude,  format = 'png', cmap='gray')
ftpng = Image.open(dir+filename2+ext3)
ftpng.save(dir+filename2+ext2)
bmpFT = Image.open(dir+filename2+ext2)
bmpFT = bmpFT.convert('L')
bmpFT.save(dir+filename2+ext2)


plt.imshow(magnitude, cmap='gray')
plt.show()

plt.imshow(phase, cmap = 'gray')
plt.show()

ift = np.fft.ifftshift(ftShift)
ift = np.fft.ifft2(ift)

plt.imshow(np.abs(ift), cmap = 'gray')
plt.show()
