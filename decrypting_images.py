#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  3 13:16:50 2021

@author: anirban
"""

import os 
import numpy as np 
import matplotlib.pyplot as plt
import nibabel as nib
from nibabel.testing import data_path
import zlib 
import dippykit as dp
#import archiver
import sys
#print(zlib.__version__)


def show_slices(slices):
   """ Function to display row of image slices """
   fig, axes = plt.subplots(1, len(slices))
   for i, slice in enumerate(slices):
       axes[i].imshow(slice.T, cmap="gray", origin="lower")
       

nii_img1 = nib.load('splitted/sub_vol0.nii.gz')
nii_img2 = nib.load('splitted/sub_vol1.nii.gz')
nii_img3 = nib.load('splitted/sub_vol2.nii.gz')
#nii_img4 = nib.load('splitted/sub_vol3.nii.gz')
#nii_img5 = nib.load('splitted/sub_vol4.nii.gz')
#nii_img6 = nib.load('splitted/sub_vol5.nii.gz')

nii_img1_data = nii_img1.get_fdata()
nii_img2_data = nii_img2.get_fdata()
nii_img3_data = nii_img3.get_fdata()
#nii_img4_data = nii_img4.get_fdata()
#nii_img5_data = nii_img5.get_fdata()
#nii_img6_data = nii_img6.get_fdata()

#maximum = np.amax()

print(np.dtype(nii_img1_data[0][0][0]))

#max_channels = np.amax([np.amax(nii_img1_data[:,:,0]), np.amax(nii_img1_data[:,:,1]), np.amax(nii_img1_data[:,:,2])])

#print (max_channels)

nii_img1_data_int = nii_img1_data.astype(np.uint16)
nii_img2_data_int = nii_img2_data.astype(np.uint16)
nii_img3_data_int = nii_img3_data.astype(np.uint16)
#nii_img4_data_int = nii_img4_data#.astype(np.uint16)
#nii_img5_data_int = nii_img5_data#.astype(np.uint16)
#nii_img6_data_int = nii_img6_data#.astype(np.uint16)

diffImage21 = nii_img1_data_int - nii_img2_data_int
diffImage32 = nii_img2_data_int - nii_img3_data_int
#diffImage43 = nii_img4_data_int - nii_img3_data_int
#diffImage54 = nii_img5_data_int - nii_img4_data_int
#diffImage65 = nii_img6_data_int - nii_img5_data_int

newImage = nii_img1_data_int + diffImage21 + diffImage32 #+ \
            #diffImage43 + diffImage54 + diffImage65

            
print(np.dtype(newImage[0][0][0]))
     
'''       
archive = archiver ( 'zip', {
        zlib : {9}
        })
    '''

#newImage = newImage.astype(np.uint8)
            
#print ("Here is the sum before : ", np.sum(newImage - nii_img6_data_int ))

newImage = newImage.copy(order='C')

#Residual_newImage = Residual_newImage.copy(order='C')

#print(nii_img1_data_C.shape)

compressImageObject = zlib.compressobj(wbits=15, memLevel=9, strategy=zlib.Z_HUFFMAN_ONLY)
blockSize = 1024

nii_img1_data_int = nii_img1_data_int.copy(order='C')
diffImage21 = diffImage21.copy(order='C')
diffImage32 = diffImage32.copy(order='C')
#diffImage43 = diffImage43.copy(order='C')
#diffImage54 = diffImage54.copy(order='C')
#diffImage65 = diffImage65.copy(order='C')

print("diffImage shape: ", diffImage21.shape)

print("Original Size : " , sys.getsizeof(nii_img1_data+nii_img2_data+\
                  nii_img3_data )) #+nii_img4_data_int+nii_img5_data_int+nii_img6_data_int))

print("Difference Original Size : " , sys.getsizeof(nii_img1_data_int+diffImage21+\
                  diffImage32 )) #+diffImage43+diffImage54+diffImage65))


maximum = 0

max1 = np.amax(nii_img1_data_int, where=~np.isnan(nii_img1_data_int), initial=0.0)
max21 = np.amax(nii_img2_data_int, where=~np.isnan(diffImage21), initial=0.0)
max32 = np.amax(nii_img3_data_int, where=~np.isnan(diffImage32), initial=0.0)



maximum = np.amax([max1, max21, max32])

print("maximum ", maximum)

print ("max21 ", max21, " max 32  ", max32)

#print (nii_img1_data_int[0][0][1])

#nii_img1_data_int[nii_img1_data_int<0] += maximum
diffImage21[diffImage21<0] += max21 
diffImage32[diffImage32<0] += max32 

#print (nii_img1_data_int[0][0][1])

compressImage = compressImageObject.compress(nii_img1_data_int) + compressImageObject.compress(diffImage21) + \
                        compressImageObject.compress(diffImage32) ##+ compressImageObject.compress(diffImage43) +\
                        #compressImageObject.compress(diffImage54) + compressImageObject.compress(diffImage65)

#compressImage_new = compressImageObject.compress(newImage)                        
                        
compressImage += compressImageObject.flush()
#compressImage_new += compressImageObject.flush()

#print(""newImage.size)

#print(nii_img1_data.size+nii_img2_data.size+nii_img3_data.size+ \
 #     nii_img4_data.size+nii_img5_data.size+nii_img6_data.size)
#compressRatio = (float(len(nii_img1_data+nii_img2_data+nii_img3_data+\
#                           nii_img4_data+nii_img5_data+nii_img6_data))-\
#                float(len(compressImage)))/float(len(nii_img1_data+nii_img2_data+nii_img3_data+\
#                           nii_img4_data+nii_img5_data+nii_img6_data))
 
print("compressed Image Size : ", float(len(compressImage)))
#print("compressed Image Size : ", float(len(compressImage_new)))
#compressRatio = (float(nii_img1_data.size+nii_img2_data.size+nii_img3_data.size+ \
 #     nii_img4_data.size+nii_img5_data.size+nii_img6_data.size)) / float(len(compressImage))

decompressedImageObject = zlib.decompressobj(wbits=+15, )
#my_file = open('compressed.dat', 'rb').read()         
#buf = my_file.read(blockSize)

image1 = np.empty((90,104,72))
image2 = np.empty((90,104,72))


#while buf:
#while compressImage:
decompressedImage = decompressedImageObject.decompress(compressImage)
decompressedImage += decompressedImageObject.flush()
decompressedImageObject = np.frombuffer(decompressedImage, dtype=np.uint16)
decompressedImageObject = np.reshape(decompressedImageObject, newshape=(-1,64,36))
splitImage = np.array_split(decompressedImageObject, 3, axis = 0)

print ("shape of split Image", splitImage[0].shape)

print ("Image 1 ", np.sum(splitImage[0] - nii_img1_data_int))
print ("Image 2 ", np.sum(splitImage[0]+splitImage[1] - nii_img2_data_int))
print ("Image 3 ", np.sum(splitImage[0]+splitImage[1]+splitImage[2] - nii_img3_data_int))

#print((decompressedImage.shape))

#print("Here is the sum of the decompressed image and the 6th image : ", np.sum(decompressedImage-nii_img6_data_int))
              
#print('Compressed: %d%%' % (100* compressRatio))

#print(new_image)

#for i in range (0, 5):

#newImage_print = newImage.astype(np.float64)

#print (np.dtype(newImage_print[0][0][0]))


#dippykit Huffman coding
directResImgUint16Max = nii_img1_data_int + diffImage21 + diffImage32

directResImgUint16Max_flatten = directResImgUint16Max.flatten()

im_encodedDirResi, stream_lengthDirResi, symbol_code_dictDirResi, symbol_prob_dictDirResi = dp.huffman_encode(directResImgUint16Max_flatten)
bitLen = 0
maxBitLen = 0
for k in symbol_code_dictDirResi:
    if maxBitLen < len(symbol_code_dictDirResi[k]):
        maxBitLen = len(symbol_code_dictDirResi[k])
    bitLen += len(symbol_code_dictDirResi[k])
totHuffmanLen = (len(symbol_code_dictDirResi)*16 + stream_lengthDirResi + len(symbol_code_dictDirResi)*maxBitLen)//8
print('Huffman Dippykit (Direct Residual):', totHuffmanLen)
print(' ')
#Making Huffman dictionary more compressible
huffmanDictKeysDirResi = symbol_code_dictDirResi.keys()
huffmanDictKeysNpArrDirResi = np.zeros(len(huffmanDictKeysDirResi), dtype='uint16')
i = 0
for k in huffmanDictKeysDirResi:
    huffmanDictKeysNpArrDirResi[i] = k
    i = i+1
# zlib Deflate
compZlibHuffmanTreeKeysDirResi = zlib.compress(huffmanDictKeysNpArrDirResi,level=9)
print('Deflate Level 9 (Huffman Tree Keys):', len(compZlibHuffmanTreeKeysDirResi))

huffmanTreeValsDirResi = symbol_code_dictDirResi.values()
huffmanTreeValListDirResi = []
for val in huffmanTreeValsDirResi:
    huffmanTreeValListDirResi.extend(val)
    huffmanTreeValListDirResi.append(2)
huffmanTreeValsNpArrDirResi = np.zeros(len(huffmanTreeValListDirResi), dtype='uint8')
i = 0
for k in huffmanTreeValListDirResi:
    if k == 0:
        huffmanTreeValsNpArrDirResi[i] = 0
    elif k == 1:
        huffmanTreeValsNpArrDirResi[i] = 1
    else:
        huffmanTreeValsNpArrDirResi[i] = 2
    i = i+1
#print(huffmanTreeValsNpArrDirResi)
    
# zlib Deflate
compHuffmanTreeValsNpArrDirResi = zlib.compress(huffmanTreeValsNpArrDirResi,level=9)
print('Deflate Level 9 (Huffman Tree Values):', len(compHuffmanTreeValsNpArrDirResi))

totHuffmanLen = len(compZlibHuffmanTreeKeysDirResi) + len(compHuffmanTreeValsNpArrDirResi) + len(im_encodedDirResi)
print('Modified Huffman Dippykit (Image):', totHuffmanLen)



'''
slice_0 = newImage_print[45, :, :]
slice_1 = newImage_print[:, 30, :]
slice_2 = newImage_print[:, :, 16]
show_slices([slice_0, slice_1, slice_2])
plt.suptitle("Center slices for EPI image") 

slice_0 = nii_img6_data[45, :, :]
slice_1 = nii_img6_data[:, 30, :]
slice_2 = nii_img6_data[:, :, 16]
show_slices([slice_0, slice_1, slice_2])
plt.suptitle("Center slices for EPI image")  

slice_0 = decompressedImage[45, :, :]
slice_1 = decompressedImage[:, 30, :]
slice_2 = decompressedImage[:, :, 16]
show_slices([slice_0, slice_1, slice_2])
plt.suptitle("Center slices for EPI image")
'''

#print(nii_img_data.header)
   


#for i in range(0, nii_img.shape[3]):
 #   print(nii_img[:,:,:,i])


