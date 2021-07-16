import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import cv2 as cv
from datetime import datetime
from scipy import ndimage, signal


def derivative_xx(image):
    width, height = image.shape[0], image.shape[1]
    len = width * height
    dxx = np.identity(width)
    for d in range(0, width - 1):
        dlen = d * width * height;
        slen = width * height;

        # print ( "dlen ", dlen )

        # print ("slen  : " ,slen)

        # print("slen - width : " , slen - width)

        for i in range(1, height - 1):
            # print (" i  ", i, "  i+1+dlen  ", i+1+dlen, " i+dlen " , i+dlen, "  i-1+dlen  ", i-1+dlen )
            dxx[i] = (image[i - 1] - 2 * image[i] + image[i + 1])

        for k in range(0, height):
            print(" k * width+dlen  ", k * width+dlen, "  k * width+1+dlen ", \
                  k * width+1+dlen, " k * width+dlen ", k * width+dlen)
            dxx[k * width+dlen] = (image[k*width+1+dlen]-image[k * width+dlen])

            print("(k+1) * width-1+dlen  ", (k+1) * width-1+dlen, " (k+1) * width-2+dlen ", \
                  (k+1) * width-2+dlen, "(k+1) * width-1+dlen ", (k+1) * width-1+dlen)
            dxx[(k+1) * width-1+dlen] = (image[(k+1) * width-2+dlen]-image[(k+1) * width-1+dlen])

            print("-----------------------------------------------")
    return dxx


def derivative_yy(image):
    width, height = image.shape[0], image.shape[1]
    len = width * height
    dyy = np.identity(width)
    for d in range(0, height - 1):
        dlen = d * width * height;
        slen = width * height;

        for i in range(1, width - 1):
            dyy[i] = (image[i - 1] - 2 * image[i] + image[i + 1])

        for k in range (0,height-1):
            dyy[k * height+dlen] = (image[k * height+1+dlen]-image[k * height+dlen])
            dyy[(k+1) * height-1+dlen] = (height[(k+1) * height-2+dlen]-image[(k+1) * height-1+dlen])
    return dyy


def read_image(path):
    image = cv.imread(path, cv.IMREAD_GRAYSCALE)
    return image


def display_image(image):
    fig, axs = plt.subplots()
    axs.imshow(image, cmap='gray')
    plt.show()

    # cv.imshow('image', image)
    # cv.waitKey(0)
    # cv.destroyWindow()


def linear_homogeneous_FSI():
    currentTime = datetime.now()

    num_of_steps = 1
    threshold = 0.2
    image = read_image('pepper.png')
    mask_image = read_image('pepper_15percentage_of_random_pixels.png')
    pixel_locations = np.loadtxt("pixel_locations", delimiter="\n", dtype=float).astype(int)

    print(pixel_locations)
    timeStepSize = 0.5
    finalImage = image - mask_image
    width, height = finalImage.shape[0], finalImage.shape[1]
    print(finalImage.shape)
    l2_norm = np.linalg.norm(image - finalImage)
    print("L2 norm = ", l2_norm)
    presentImage = finalImage
    prevImage = np.identity(width)
    currentImage = finalImage

    #while l2_norm > dxx:
    for n in range(0, num_of_steps):
        dxx = ndimage.correlate1d(presentImage, [1, -2, 1], axis=0)  # derivative_xx(finalImage)
        dyy = ndimage.correlate1d(currentImage, [1, -2, 1], axis=1)  # derivative_yy(finalImage)
        # print("calculated dxx and dyy")
        #display_image(dxx)
        #display_image(dyy)
        #dImage = np.hypot(dxx, dyy)
        dImage = (dxx+dyy)
        #print("Dimages  = ", ((dImageOld - dImage) ** 2).mean())
        #display_image(dImage)
        alpha = (4 * n + 2) / (2 * n + 3)
        currentImage = presentImage
        traceIndex = 0
#        for i in range (0 ,width):
#         if n == 0:
#             prevImage = presentImage
        presentImage = alpha * (presentImage + (timeStepSize * (dImage))) + (1 - alpha) * prevImage
        prevImage = currentImage

    l2_norm = np.linalg.norm(image - presentImage)
    print("L2 norm = ", l2_norm)

    mse = ((image - presentImage) ** 2).mean()
    print("Mean Squared Error  = ", mse)

    print("Run Time = ", datetime.now() - currentTime)
    display_image(presentImage)
    return 0

def EED_FSI():
    num_of_steps = 5
    currentTime = datetime.now()
    image = read_image('pepper.png')
    mask_image = read_image('pepper_15percentage_of_random_pixels.png')
    pixel_locations = np.loadtxt("pixel_locations", delimiter="\n", dtype=float).astype(int)
    finalImage = image - mask_image
    presentImage = finalImage

    t = np.linspace(-10,10,30)
    bump = np.exp(-0.1 * t ** 2)
    bump /= np.trapz(bump)  # normalize the integral to 1
    # make a 2-D kernel out of it
    kernel = bump[:, np.newaxis] * bump[np.newaxis, :]

    # u_sigma = signal.fftconvolve(finalImage, kernel[:,:,np.newaxis], mode='same')
    width, height = u_sigma.shape[0], u_sigma.shape[1]

    presentImage = u_sigma
    prevImage = np.identity(width)
    currentImage = u_sigma
    for n in range(0, num_of_steps):
        u_sigma = signal.fftconvolve(finalImage, kernel[:, :, np.newaxis], mode='same')
        dxx_gaussian = ndimage.correlate1d(u_sigma, [1, -2, 1], axis=0)  # derivative_xx(finalImage)
        dyy_gaussian = ndimage.correlate1d(u_sigma, [1, -2, 1], axis=1)  # derivative_yy(finalImage)
        dxx = ndimage.correlate1d(finalImage, [1, -2, 1], axis=0)  # derivative_xx(finalImage)
        dyy = ndimage.correlate1d(finalImage, [1, -2, 1], axis=1)  # derivative_yy(finalImage)
        # print("calculated dxx and dyy")
        # display_image(dxx)
        # display_image(dyy)
        # dImage = np.hypot(dxx, dyy)
        dImage = (dxx + dyy)
        # print("Dimages  = ", ((dImageOld - dImage) ** 2).mean())
        # display_image(dImage)
        alpha = (4 * n + 2) / (2 * n + 3)
        currentImage = presentImage
        traceIndex = 0
        #        for i in range (0 ,width):
        #         if n == 0:
        #             prevImage = presentImage
        presentImage = alpha * (presentImage + (timeStepSize * (dImage))) + (1 - alpha) * prevImage
        prevImage = currentImage



    print("Run Time = ", datetime.now() - currentTime)
    display_image(presentImage)
    return 0