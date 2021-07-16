import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import cv2 as cv
from datetime import datetime
from scipy import ndimage, signal


def derivative_xx(image, width, height, stepsize):
    #width, height = image.shape[0], image.shape[1]
    len = width * height
    dxx = np.ones(len)


    # out_array
    for i in range(1, len - 1):
        # print (" i  ", i, "  i+1+dlen  ", i+1+dlen, " i+dlen " , i+dlen, "  i-1+dlen  ", i-1+dlen )
        dxx[i] = (image[i + 1] - 2 * image[i] + image[i - 1])/(stepsize ** 2)

    for k in range(0, height):
        # print(" k * width+dlen  ", k * width+dlen, "  k * width+1+dlen ", \
        #             k * width+1+dlen, " k * width+dlen ", k * width+dlen)
        dxx[k * width] = (image[k * width + 1] - image[k * width])

        # print("(k+1) * width-1+dlen  ", (k+1) * width-1+dlen, " (k+1) * width-2+dlen ", \
        #           (k+1) * width-2+dlen, "(k+1) * width-1+dlen ", (k+1) * width-1+dlen)
        dxx[(k + 1) * width - 1] = (image[(k + 1) * width - 2] - image[(k + 1) * width - 1])

        #print("-----------------------------------------------")
    return dxx


def derivative_yy(image, width, height, stepsize):
    #width, height = image.shape[0], image.shape[1]
    len = width * height
    dyy = np.ones(len)
    for i in range(width, len - width):
        dyy[i] = (image[i + width] - 2 * image[i] + image[i - width])/(stepsize ** 2)

    for k in range(0, width):
        dyy[k] = (image[k + width] - image[k])
        dyy[(height - 1) * width + k ] = (image[(height - 2) * width + k] - image[k + width * (height - 1)])
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

    num_of_steps = 100
    stepsize = 1
    threshold = 0.2
    image = read_image('pepper.png')
    mask_image = read_image('test_mask/pepper_mask_15percentage_of_random_pixels.png')
    pixel_locations = np.loadtxt('test_mask/pixel_locations', delimiter="\n", dtype=float).astype(int)

    print(pixel_locations)
    timeStepSize = 0.25
    #finalImage = image - mask_image
    #width, height = finalImage.shape[0], finalImage.shape[1]
    #print(finalImage.shape)

    #converting to 1d
    original_width, original_height = image.shape[0], image.shape[1]
    masked_width, masked_height =  mask_image.shape[0], mask_image.shape[1]

    image_1d = image.reshape(original_width*original_height)
    maskedImage_1d = mask_image.reshape(masked_width*masked_height)

    l2_norm = np.linalg.norm(image - mask_image)
    print("L2 norm = ", l2_norm)

    presentImage = maskedImage_1d
    prevImage = maskedImage_1d
    currentImage = maskedImage_1d

    print("masked image  =  ", presentImage[7550], "   Original Image =  ", image_1d[7550])

    # while l2_norm > threshold:
    for n in range(0, num_of_steps):
        dxx = derivative_xx(maskedImage_1d, masked_width, masked_height, stepsize)  # ndimage.correlate1d(presentImage, [1, -2, 1], axis=0)
        dyy = derivative_yy(maskedImage_1d, masked_width, masked_height, stepsize)  # ndimage.correlate1d(currentImage, [1, -2, 1], axis=1)

        print("dxxleft image  =  ", dxx[7550-1])
        print("dxxtop image  =  ", dxx[7550-masked_width])
        print("dxx image  =  ", dxx[7550])
        print("dxxright image  =  ", dxx[7550 + 1])
        print("dxxbottom image  =  ", dxx[7550 + masked_width])

        print("dyy image  =  ", dxx[7550])

        # print("calculated dxx and dyy")
        # display_image(dxx)
        # display_image(dyy)
        # dImage = np.hypot(dxx, dyy)
        dImage = (dxx + dyy)
        # print("Dimages  = ", ((dImageOld - dImage) ** 2).mean())
        # display_image(dImage)
        alpha = (4 * n + 2) / (2 * n + 3)
        #print(alpha)
        traceIndex = 0
        #print(masked_width*masked_height)
        for i in range (0, masked_width*masked_height-1):
            if n == 0:
                prevImage[i] = presentImage[i]

            currentImage[i] = presentImage[i]

            if (i == pixel_locations[traceIndex]):
                #print("here at ", traceIndex)
                traceIndex+=1
                #print(traceIndex)
                continue

            presentImage[i] = alpha * (presentImage[i] + (timeStepSize * (dImage[i]))) + (1 - alpha) * prevImage[i]
            prevImage[i] = currentImage[i]

    print("masked image  =  ", presentImage[7550], "   Original Image =  ", image_1d[7550])

    l2_norm = np.linalg.norm(image_1d - presentImage)
    print("L2 norm = ", l2_norm)

    mse = ((image_1d - presentImage) ** 2).mean()
    print("Mean Squared Error  = ", mse)

    print_image = presentImage.reshape(original_width,original_height)

    print("Run Time = ", datetime.now() - currentTime)
    display_image(print_image)
    return 0


def EED_FSI():
    num_of_steps = 5
    currentTime = datetime.now()
    image = read_image('pepper.png')
    mask_image = read_image('pepper_15percentage_of_random_pixels.png')
    pixel_locations = np.loadtxt("pixel_locations", delimiter="\n", dtype=float).astype(int)
    finalImage = image - mask_image
    presentImage = finalImage

    t = np.linspace(-10, 10, 30)
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
