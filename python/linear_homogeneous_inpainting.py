import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import cv2 as cv


def read_image(path):
    image = cv.imread(path, cv.IMREAD_GRAYSCALE )
    return image

def display_image(image):
    cv.imshow('image', image)
    cv.waitKey(0)
    cv.destroyWindow()


def linear_homogeneous_FSI():

    threshold = 0.2
    image = read_image('pepper.png')
    mask_image = read_image('pepper_15percentage_of_random_pixels.png')

    finalImage = image - mask_image

    print(finalImage.shape)

    l2_norm = np.linalg.norm(image-finalImage)
    print("L2 norm = " , l2_norm)

    while l2_norm>threshold:
        

    #display_image(finalImage)
    return 0
