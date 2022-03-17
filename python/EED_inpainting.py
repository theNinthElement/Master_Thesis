import numpy as np
import matplotlib.pyplot as plt
import cv2 as cv
from datetime import datetime
from scipy import ndimage, signal


def derivative_xx(image, width, height, stepsize):
    # width, height = image.shape[0], image.shape[1]
    len = width * height
    dxx = np.zeros(len)
    # out_array
    for i in range(1, len - 1):
        # print (" i  ", i, "  i+1+dlen  ", i+1+dlen, " i+dlen " , i+dlen, "  i-1+dlen  ", i-1+dlen )
        dxx[i] = (image[i + 1] - 2 * image[i] + image[i - 1]) / (stepsize * stepsize)

    for k in range(0, height):
        # print(" k * width+dlen  ", k * width+dlen, "  k * width+1+dlen ", \
        #             k * width+1+dlen, " k * width+dlen ", k * width+dlen)
        dxx[k * width] = (image[k * width + 1] - image[k * width]) / (stepsize * stepsize)

        # print("(k+1) * width-1+dlen  ", (k+1) * width-1+dlen, " (k+1) * width-2+dlen ", \
        #           (k+1) * width-2+dlen, "(k+1) * width-1+dlen ", (k+1) * width-1+dlen)
        dxx[(k + 1) * width - 1] = (image[(k + 1) * width - 2] - image[(k + 1) * width - 1]) / (stepsize * stepsize)

        # print("-----------------------------------------------")
    return dxx


def derivative_yy(image, width, height, stepsize):
    # width, height = image.shape[0], image.shape[1]
    len = width * height
    dyy = np.zeros(len)
    for i in range(width, len - width):
        dyy[i] = (image[i + width] - 2 * image[i] + image[i - width]) / (stepsize * stepsize)

    for k in range(0, width):
        dyy[k] = (image[k + width] - image[k]) / (stepsize * stepsize)
        dyy[(height - 1) * width + k] = (image[(height - 2) * width + k] - image[(height - 1) * width + k]) / (
                    stepsize * stepsize)
    return dyy


def dervForwX(image, width, height, stepsize):
    len = width * height
    outAr = np.zeros(len)
    for i in range(0, len):
        if (i + 1) % width == 0:
            outAr[i] = 0
        else:
            outAr[i] = (image[i + 1] - image[i]) / stepsize
    return outAr


def dervForwY(image, width, height, stepsize):
    len = width * height
    outAr = np.zeros(len)
    for i in range(0, len-width):
        outAr[i] = (image[i+width]-image[i])/stepsize;
    for i in range(width*(height-1), len):
        outAr[i] = 0
    return outAr


def dervBackX(image, width, height, stepsize):
    len = width * height
    outAr = np.zeros(len)
    for i in range(0, len):
        if (i + 1) % width == 0:
            outAr[i] = 0
        else:
            outAr[i] = (image[i] - image[i-1]) / stepsize
    return outAr


def dervBackY(image, width, height, stepsize):
    len = width * height
    outAr = np.zeros(len)
    for i in range(0, width):
       outAr[i] = 0
    for i in range(width, len):
        outAr[i] = (image[i]-image[i-width])/stepsize;
    return outAr


# def Convolute_X (image, width, height, kernel):
#
#     for i in range (0, height):
#         offset = 0
#         for j in range (0)
#
#     return convolution_X


def read_image(path):
    image = cv.imread(path, cv.IMREAD_GRAYSCALE)
    return image


def display_image(image):
    fig, axs = plt.subplots()
    axs.imshow(image, cmap='gray')
    plt.show()


def EED_FSI():
    num_of_steps = 5
    currentTime = datetime.now()
    image = read_image('pepper.png')
    mask_image = read_image('test_mask/pepper_mask_15percentage_of_random_pixels.png')
    pixel_locations = np.loadtxt("test_mask/pixel_locations", delimiter="\n", dtype=float).astype(int)
    finalImage = image - mask_image
    presentImage = finalImage

    # converting to 1d
    original_width, original_height = image.shape[0], image.shape[1]
    masked_width, masked_height = mask_image.shape[0], mask_image.shape[1]

    image_1d = image.reshape(original_width * original_height)
    image_1d = image_1d.astype(float)
    maskedImage_1d = mask_image.reshape(masked_width * masked_height)
    maskedImage_1d = maskedImage_1d.astype(float)

    presentImage = maskedImage_1d

    # t = np.linspace(-10, 10, 30)
    # bump = np.exp(-0.1 * t ** 2)
    # bump /= np.trapz(bump)  # normalize the integral to 1
    # make a 1-D kernel out of it
    kernel = ndimage.gaussian_filter1d(np.float_([0, 1, 0]), 1)  # bump[:, np.newaxis] * bump[np.newaxis, :]

    # u_sigma = signal.fftconvolve(finalImage, kernel[:,:,np.newaxis], mode='same')
    # width, height = u_sigma.shape[0], u_sigma.shape[1]


    prevImage = maskedImage_1d
    currentImage = maskedImage_1d

    d11 = np.zeros(masked_width * masked_height)
    d12 = np.zeros(masked_width * masked_height)
    d21 = np.zeros(masked_width * masked_height)
    d22 = np.zeros(masked_width * masked_height)

    for n in range(0, num_of_steps):
        # u_sigma = signal.fftconvolve(finalImage, kernel[:, :, np.newaxis], mode='same')

        convolove_x = ndimage.gaussian_filter1d(presentImage, 1)

        dxx_gaussian = derivative_xx(convolove_x, masked_width, masked_height, 1)
        dyy_gaussian = derivative_yy(convolove_x, masked_width, masked_height, 1)
        dxx = derivative_xx(presentImage, masked_width, masked_height, 1)
        dyy = derivative_yy(presentImage, masked_width, masked_height, 1)

        dervForX = dervForwX(presentImage, masked_width, masked_height, 1);
        dervForY = dervForwY(presentImage, masked_width, masked_height, 1);
        dervBacX = dervBackX(presentImage, masked_width, masked_height, 1);
        dervBacY = dervBackY(presentImage, masked_width, masked_height, 1);

        for i in range (0, masked_width * masked_height):
            if n == 0:
                prevImage[i] = presentImage[i];
            currentImage[i] = presentImage[i]

            norm_i_square = dxx_gaussian[i] * dxx_gaussian[i] + dyy_gaussian[i] * dyy_gaussian[i];
            normi = np.sqrt(norm_i_square)
            d1 = 0.0
            d2 = 0.0
            v1 = [0.0,0.0]
            v2 = [0.0,0.0]

            if normi == 0:
                v1[0] = 1;
                v1[1] = 0;

                v2[0] = 0;
                v2[1] = 1;
            else:
                v1[0] = dxx_gaussian[i] / normi
                v1[1] = dyy_gaussian[i] / normi

                v2[0] = dyy_gaussian[i] / normi
                v2[1] = -dxx_gaussian[i] / normi

            d1 = 1.0/np.sqrt(1 + norm_i_square/(0.5**2));
            d2 = 1;

            d11[i] = d1 * v1[0] * v1[0] + d2 * v2[0] * v2[0];
            d12[i] = (d1 * v1[0] * v1[1] + d2 * v2[0] * v2[1]) * dyy[i];
            d21[i] = (d1 * v1[1] * v1[0] + d2 * v2[1] * v2[0]) * dxx[i];
            d22[i] = d1 * v1[1] * v1[1] + d2 * v2[1] * v2[1];

        sumForX = sumForwX(d11, imgWidth, imgHeight);
        sumForY = sumForwY(d22, imgWidth, imgHeight);
        sumBacX = sumBackX(d11, imgWidth, imgHeight);
        sumBacY = sumBackY(d22, imgWidth, imgHeight);


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
