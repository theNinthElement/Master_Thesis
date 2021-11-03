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
import gzip

#print(zlib.__version__)


def show_slices(slices):
   """ Function to display row of image slices """
   fig, axes = plt.subplots(1, len(slices))
   for i, slice in enumerate(slices):
       axes[i].imshow(slice.T, cmap="gray", origin="lower")

nii_img1 = nib.load('fMRI_splitted/sub_vol0.nii.gz')
nii_img2 = nib.load('fMRI_splitted/sub_vol1.nii.gz')
nii_img3 = nib.load('fMRI_splitted/sub_vol2.nii.gz')
nii_img4 = nib.load('fMRI_splitted/sub_vol3.nii.gz')
nii_img5 = nib.load('fMRI_splitted/sub_vol4.nii.gz')
nii_img6 = nib.load('fMRI_splitted/sub_vol5.nii.gz')
#nii_img  = nib.load('someones_anatomy.nii.gz')
#nii_img  = nib.load('someones_epi.nii.gz')

#nii_data = nii_img.get_fdata()

#print (len(nii_data))

#nii_img_data = nii_img.get_fdata()
#epi_img_data.shape
#print(nii_img_data.shape)


nii_img1_data = nii_img1.get_fdata()
nii_img2_data = nii_img2.get_fdata()
nii_img3_data = nii_img3.get_fdata()
nii_img4_data = nii_img4.get_fdata()
nii_img5_data = nii_img5.get_fdata()
nii_img6_data = nii_img6.get_fdata()

print(np.dtype(nii_img1_data[0][0][0]))

nii_img1_data_int = nii_img1_data#.astype(np.uint8)
nii_img2_data_int = nii_img2_data#.astype(np.uint8)
nii_img3_data_int = nii_img3_data#.astype(np.uint8)
nii_img4_data_int = nii_img4_data#.astype(np.uint8)
nii_img5_data_int = nii_img5_data#.astype(np.uint8)
nii_img6_data_int = nii_img6_data#.astype(np.uint8)

diffImage21 = nii_img2_data_int - nii_img1_data_int
diffImage32 = nii_img3_data_int - nii_img2_data_int
diffImage43 = nii_img4_data_int - nii_img3_data_int
diffImage54 = nii_img5_data_int - nii_img4_data_int
diffImage65 = nii_img6_data_int - nii_img5_data_int

newImage = nii_img1_data_int + diffImage21 + diffImage32 + \
            diffImage43 + diffImage54 + diffImage65
            
print(np.dtype(newImage[0][0][0]))
            
#newImage = newImage.astype(np.uint8)
            
print ("Here is the sum before : ", np.sum(newImage - nii_img6_data_int ))

newImage = newImage.copy(order='C')
nii_img1_data_C = nii_img1_data_int.copy(order='C')

#print(nii_img1_data_C.shape)

compressImageObject = zlib.compressobj(strategy=zlib.Z_HUFFMAN_ONLY, wbits=15, memLevel=9)
blockSize = 1024
compressImage = compressImageObject.compress(newImage)
compressImage += compressImageObject.flush()

print(newImage.size)

#print(nii_img1_data.size+nii_img2_data.size+nii_img3_data.size+ \
 #     nii_img4_data.size+nii_img5_data.size+nii_img6_data.size)
#compressRatio = (float(len(nii_img1_data+nii_img2_data+nii_img3_data+\
#                           nii_img4_data+nii_img5_data+nii_img6_data))-\
#                float(len(compressImage)))/float(len(nii_img1_data+nii_img2_data+nii_img3_data+\
#                           nii_img4_data+nii_img5_data+nii_img6_data))
 
print(float(len(compressImage)))
#compressRatio = (float(nii_img1_data.size+nii_img2_data.size+nii_img3_data.size+ \
 #     nii_img4_data.size+nii_img5_data.size+nii_img6_data.size)) / float(len(compressImage))

decompressedImageObject = zlib.decompressobj(wbits=+15, )
#my_file = open('compressed.dat', 'rb').read()         
#buf = my_file.read(blockSize)

#while buf:
decompressedImage = decompressedImageObject.decompress(compressImage)
#    buf = my_file.read(blockSize)
    
decompressedImage += decompressedImageObject.flush()
print(float(len(decompressedImage)))
decompressedImage = np.frombuffer(decompressedImage, dtype=np.float64)
decompressedImage = np.reshape(decompressedImage, newshape=(90,104,72))
#decompressedImage = decompressedImage.astype(np.float64)
#decompressedImage = decompressedImage.decode('utf-8')

print((decompressedImage.shape))

print("Here is the sum of the decompressed image and the 6th image : ", np.sum(decompressedImage-nii_img6_data_int))
              
#print('Compressed: %d%%' % (100* compressRatio))

#print(new_image)

#for i in range (0, 5):

newImage_print = newImage.astype(np.float64)

print (np.dtype(newImage_print[0][0][0]))



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

#print(nii_img_data.header)
   


#for i in range(0, nii_img.shape[3]):
 #   print(nii_img[:,:,:,i])


