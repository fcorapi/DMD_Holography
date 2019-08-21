#gaussianFit2D
#This program fits a gaussian to an image. Could be useful.
#Adapted from Scipy Cookbook (https://scipy-cookbook.readthedocs.io/items/FittingData.html)
#Frank Corapi (fcorapi@uwaterloo.ca)
#June 27th, 2019

import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.image as img
from scipy import optimize

def gaussian(a, x_0, y_0, wx, wy):
    wx = float(wx)
    wy = float(wy)
    return lambda x,y: a*np.exp(-np.power(((x-x_0)/wx),2)-np.power(((y-y_0)/wy),2))

def moments(data):
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X * data).sum() / total
    y = (Y * data).sum() / total
    col = data[:, int(y)]
    width_x = np.sqrt(np.abs((np.arange(col.size) - y) ** 2 * col).sum() / col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(np.abs((np.arange(row.size) - x) ** 2 * row).sum() / row.sum())
    height = data.max()
    return height, x, y, width_x, width_y

def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) -
                                 data)
    p, success = optimize.leastsq(errorfunction, params)
    return p

dir = 'C:\Users\Franky\Desktop\UofT Summer 2019\ImagePlane\\'
filename = 'Gaussian2_Cropped'
ext1 = '.bmp'

initialImg = Image.open(dir+filename+ext1)
plt.imshow(initialImg)
plt.show()