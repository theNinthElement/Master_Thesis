import numpy as np
import matplotlib.pyplot as plt
import cv2 as cv
from datetime import datetime
from scipy import ndimage, signal
import nibabel as nib


def charbonnier_diff(s_square, contrastParam):
    g = 1.0 / np.sqrt(1 + s_square / (contrastParam ** 2))
    return g


def gaussian1D_kernel(gkernel, kernelLength, sigma):
    sum = 0.0
    s = 2.0 * (sigma ** 2)

    # for i in range(0, kernelLength):
    #     gkernel

    for x in range(-(kernelLength - 1) / 2, (kernelLength - 1) / 2):
        r = x ** 2
        gkernel[x + (kernelLength - 1) / 2] = np.exp(-r / s) / np.sqrt(np.pi * s)
        sum += gkernel[x + (kernelLength - 1) / 2]

    for i in range(0, kernelLength):
        gkernel[i] /= sum

    return gkernel


def derivative3DX(image, width, height, depth,
                  stepsize):
    len = width * height * depth
    outputArray = np.zeros(len)

    slen = width * height

    for d in range(0, depth):
        dlen = d * width * height

        for i in range(1, slen - 1):
            outputArray[i + dlen] = (image[i + 1 + dlen] - image[i - 1 + dlen]) / (2 * (stepsize ** 2))
        for k in range(0, height):
            outputArray[k * width + dlen] = (image[k * width + 1 + dlen] - \
                                             image[k * width + dlen]) / (2 * stepsize)
            outputArray[(k + 1) * width - 1 + dlen] = (image[(k + 1) * width - 1 + dlen] - \
                                                       image[(k + 1) * width - 2 + dlen]) / (2 * (stepsize ** 2))
    return outputArray


def derivative3DY(image, width, height, depth,
                  stepsize):
    len = width * height * depth
    outputArray = np.zeros(len)

    slen = width * height

    for d in range(0, depth):
        dlen = d * width * height

        for i in range(1, slen - 1):
            outputArray[i + dlen] = (image[i + width + dlen] - image[i - width + dlen]) / (2 * (stepsize ** 2))
        for k in range(0, height):
            outputArray[k + dlen] = (image[k * width + dlen] - \
                                     image[k + dlen]) / (2 * stepsize)
            outputArray[(height - 1) * width + k + dlen] = (image[(height - 1) * width + k + dlen] - \
                                                            image[(height - 1) * width + k + dlen]) / (
                                                                   2 * (stepsize ** 2))
    return outputArray


def derivative3DZ(image, width, height, depth,
                  stepsize):
    len = width * height * depth
    outputArray = np.zeros(len)

    slen = width * height
    for i in range(slen, len - slen):
        outputArray[i] = (image[i + slen] - image[i - slen]) / (2 * (stepsize ** 2))
    for h in range(0, height):
        hlen = h * width
        for k in range(0, width):
            outputArray[k + hlen] = (image[k * slen + hlen] - \
                                     image[k + hlen]) / (2 * (stepsize ** 2))
            outputArray[len - slen + k + hlen] = (image[len - slen + k + hlen] - \
                                                  image[len - 2 * slen + k + hlen]) / (2 * (stepsize ** 2))
    return outputArray


def derivative3DXX(image, width, height, depth,
                   stepsize):
    len = width * height * depth
    outputArray = np.zeros(len)

    slen = width * height
    for d in range(0, depth):
        dlen = d * width * height

        for i in range(1, slen - 1):
            outputArray[i + dlen] = (image[i + 1 + dlen] - image[i - 1 + dlen]) / (stepsize ** 2)
        for k in range(0, height):
            outputArray[k * width + dlen] = (image[k * width + 1 + dlen] - \
                                             image[k * width + dlen]) / stepsize ** 2
            outputArray[(k + 1) * width - 1 + dlen] = (image[(k + 1) * width - 1 + dlen] - \
                                                       image[(k + 1) * width - 2 + dlen]) / stepsize ** 2
    return outputArray


def derivative3DYY(image, width, height, depth,
                   stepsize):
    len = width * height * depth
    outputArray = np.zeros(len)

    slen = width * height

    for d in range(0, depth):
        dlen = d * width * height

        for i in range(width, slen - width):
            outputArray[i + dlen] = (image[i + width + dlen] - 2 * image[i + dlen] + \
                                     image[i - width + dlen]) / (stepsize ** 2)
        for k in range(0, width):
            outputArray[k + dlen] = (image[k * width + dlen] - \
                                     image[k + dlen]) / (stepsize ** 2)
            outputArray[(height - 1) * width + k + dlen] = (image[(height - 2) * width + k + dlen] - \
                                                            image[(height - 1) * width + k + dlen]) / (
                                                                   2 * (stepsize ** 2))
    return outputArray


def derivative3DZZ(image, width, height, depth,
                   stepsize):
    len = width * height * depth
    outputArray = np.zeros(len)

    slen = width * height
    for i in range(slen, len - slen):
        outputArray[i] = (image[i - slen] - 2 * image[i] + image[i - slen]) / (stepsize ** 2)
    for h in range(0, height):
        hlen = h * width
        for k in range(0, width):
            outputArray[k + hlen] = (image[k * slen + hlen] - \
                                     image[k + hlen]) / (2 * (stepsize ** 2))
            outputArray[len - slen + k + hlen] = (image[len - 2 * slen + k + hlen] - \
                                                  image[len - slen + k + hlen]) / (stepsize ** 2)
    return outputArray


def convolution3DX(image, width, height, depth,
                   stepsize, gaussian_kernel, kernelLength):
    len = width, height, depth
    kCenter = kernelLength >> 1

    convolvedArr = np.zeros(len)
    endIndex = width - kCenter;

    inputArrIndex = 0
    convolvedArrIndex = 0

    for z in range(0, depth):
        for i in range(0, height):
            kOffset = 0

            for j in range(0, kCenter):
                convolvedArr[convolvedArrIndex] = 0
                m = 0
                for k in range(kCenter + kOffset, 0, -1):
                    convolvedArr[convolvedArrIndex] += \
                        image[inputArrIndex + m] * gaussian_kernel[k]
                    m += 1
                convolvedArrIndex += 1
                kOffset += 1

            for j in range(kCenter, endIndex):
                convolvedArr[convolvedArrIndex] = 0
                m = 0
                for k in range(kernelLength - 1, 0, -1):
                    convolvedArr[convolvedArrIndex] += \
                        image[inputArrIndex + m] * gaussian_kernel[k]
                    m += 1
                inputArrIndex += 1
                convolvedArrIndex += 1
            kOffset = 1

            for j in range(endIndex, width):
                convolvedArr[convolvedArrIndex] = 0
                m = 0
                for k in range(kernelLength - 1, kOffset, -1):
                    convolvedArr[convolvedArrIndex] += \
                        image[inputArrIndex + m] * gaussian_kernel[k]
                    m += 1
                inputArrIndex += 1
                convolvedArrIndex += 1
                kOffset += 1

            inputArrIndex += kCenter

    return convolvedArr


def convolution3DY(image, width, height, depth,
                   stepsize, gaussian_kernel, kernelLength):
    len = width, height, depth
    kCenter = kernelLength >> 1
    sum = np.zeros(width)
    convolvedArr = np.zeros(len)
    endIndex = height - kCenter

    for z in range(0, depth):
        inputArrIndex = 0
        currentRowNumber = 0

        kOffset = 0
        for i in range(0, kCenter):

            for k in range(kCenter + kOffset, 0, -1):
                for j in range(0, width):
                    sum[j] += image[inputArrIndex + \
                                    z * width * height] * gaussian_kernel[k]
                    inputArrIndex += 1
            for n in range(0, width):
                convolvedArr[n + i * width + \
                             z * width * height] = sum[n]
                sum[n] = 0
            inputArrIndex = currentRowNumber
            kOffset += 1

        for i in range(kCenter, endIndex):

            for k in range(kernelLength - 1, 0, -1):
                for j in range(0, width):
                    sum[j] += image[inputArrIndex + \
                                    z * width * height] * gaussian_kernel[k]
                    inputArrIndex += 1
            for n in range(0, width):
                convolvedArr[n + i * width + \
                             z * width * height] = sum[n]
                sum[n] = 0
            currentRowNumber += width
            inputArrIndex = currentRowNumber
        kOffset = 1
        for i in range(endIndex, height):

            for k in range(kernelLength - 1, kOffset, -1):
                for j in range(0, width):
                    sum[j] += image[inputArrIndex + \
                                    z * width * height] * gaussian_kernel[k]
                    inputArrIndex += 1
            for n in range(0, width):
                convolvedArr[n + i * width + \
                             z * width * height] = sum[n]
                sum[n] = 0
            currentRowNumber += width
            inputArrIndex = currentRowNumber
            kOffset += 1

    return convolvedArr


def convolution3DZ(image, width, height, depth,
                   stepsize, gaussian_kernel, kernelLength):
    len = width, height, depth
    kCenter = kernelLength >> 1
    sum = np.zeros(width)
    convolvedArr = np.zeros(len)
    endIndex = height - kCenter

    for h in range(0, height):
        inputArrIndex = 0
        currentDepthNumber = 0

        kOffset = 0
        for i in range(0, kCenter):
            for k in range(kCenter + kOffset, 0, -1):
                for j in range(0, width):
                    sum[j] += image[inputArrIndex + \
                                    h * width * height] * gaussian_kernel[k]
                    inputArrIndex += 1
                inputArrIndex = inputArrIndex - width + width * height
            for n in range(0, width):
                convolvedArr[n + i * width * height + \
                             h * width] = sum[n]
                sum[n] = 0
            inputArrIndex = currentDepthNumber
            kOffset += 1

        for i in range(kCenter, endIndex):
            for k in range(kernelLength - 1, 0, -1):
                for j in range(0, width):
                    sum[j] += image[inputArrIndex + \
                                    h * width] * gaussian_kernel[k]
                    inputArrIndex += 1
                inputArrIndex = inputArrIndex - width + width * height

            for n in range(0, width):
                convolvedArr[n + i * width * height + \
                             h * width] = sum[n]
                sum[n] = 0
            currentDepthNumber += width
            inputArrIndex = currentDepthNumber
        kOffset = 1
        for i in range(endIndex, height):

            for k in range(kernelLength - 1, kOffset, -1):
                for j in range(0, width):
                    sum[j] += image[inputArrIndex + \
                                    h * width] * gaussian_kernel[k]
                    inputArrIndex += 1
                inputArrIndex = inputArrIndex - width + width * height
            for n in range(0, width):
                convolvedArr[n + i * width * height + \
                             h * width] = sum[n]
                sum[n] = 0
            currentDepthNumber += width
            inputArrIndex = currentDepthNumber
            kOffset += 1

    return convolvedArr


def read_image(path):
    image = cv.imread(path, cv.IMREAD_GRAYSCALE)
    return image


def display_image(image):
    fig, axs = plt.subplots()
    axs.imshow(image, cmap='gray')
    plt.show()


def show_slices(slices):
    """ Function to display row of image slices """
    fig, axes = plt.subplots(1, len(slices))
    for i, slice in enumerate(slices):
        axes[i].imshow(slice.T, cmap="gray", origin="lower")


def EED_FSI():
    num_of_steps = 5
    currentTime = datetime.now()

    timeStep = 1

    randomState = np.random.RandomState(12345)

    nii_img = nib.load('./../fMRI_splitted/sub_vol2.nii.gz')

    image = nii_img.get_fdata()
    print(image.shape)
    print(image.shape[0] * image.shape[1] * image.shape[2])

    # selecting 80% random integers
    randomPixels = randomState.randint(1,
                                     (image.shape[0] * image.shape[1] * image.shape[2]),
                                     539136)
                                     #(80 * (image.shape[0] * image.shape[1] * image.shape[2]) / 100))

    #print(len(randomPixels))
    # image = read_image('pepper.png')
    # mask_image = read_image(
    #   'test_mask/pepper_mask_15percentage_of_random_pixels.png')
    # pixel_locations = np.loadtxt("test_mask/pixel_locations", delimiter="\n", dtype=float).astype(int)
    # finalImage = image - mask_image

   # presentImage = finalImage

    kernelSize = 3
    stepSize = 1
    sigma = 1.0

    # converting to 1d
    original_width, original_height, original_depth = image.shape[0], image.shape[1], image.shape[2]
    masked_width, masked_height, masked_depth = image.shape[0], image.shape[1], image.shape[2]

    mask_image = image.reshape(original_width * original_height * original_depth)

    # removing the selected pixels
    randArrTraceIndex = 0
    for i in range(0, image.shape[0] * image.shape[1] * image.shape[2] ):
        if i == randomPixels[randArrTraceIndex]:
            mask_image[i] == 0
            randArrTraceIndex +=1

    image_1d = image.reshape(original_width * original_height * original_depth)
    image_1d = image_1d.astype(float)
    maskedImage_1d = mask_image#mask_image.reshape(masked_width * masked_height * masked_depth)
    maskedImage_1d = maskedImage_1d.astype(float)

    presentImage = maskedImage_1d

    # t = np.linspace(-10, 10, 30)
    # bump = np.exp(-0.1 * t ** 2)
    # bump /= np.trapz(bump)  # normalize the integral to 1
    # make a 1-D kernel out of it
    kernel = ndimage.gaussian_filter1d(np.float_([0, 1, 0]), 1)  # bump[:, np.newaxis] * bump[np.newaxis, :]

    # u_sigma = signal.fftconvolve(finalImage, kernel[:,:,np.newaxis], mode='same')
    # width, height = u_sigma.shape[0], u_sigma.shape[1]

    d1111 = np.zeros([maskedImage_1d.shape[0]])
    d1112 = np.zeros([maskedImage_1d.shape[0]])
    d1113 = np.zeros([maskedImage_1d.shape[0]])
    d1121 = np.zeros([maskedImage_1d.shape[0]])
    d1122 = np.zeros([maskedImage_1d.shape[0]])
    d1123 = np.zeros([maskedImage_1d.shape[0]])
    d1131 = np.zeros([maskedImage_1d.shape[0]])
    d1132 = np.zeros([maskedImage_1d.shape[0]])
    d1133 = np.zeros([maskedImage_1d.shape[0]])

    d1211 = np.zeros([maskedImage_1d.shape[0]])
    d1212 = np.zeros([maskedImage_1d.shape[0]])
    d1213 = np.zeros([maskedImage_1d.shape[0]])
    d1221 = np.zeros([maskedImage_1d.shape[0]])
    d1222 = np.zeros([maskedImage_1d.shape[0]])
    d1223 = np.zeros([maskedImage_1d.shape[0]])
    d1231 = np.zeros([maskedImage_1d.shape[0]])
    d1232 = np.zeros([maskedImage_1d.shape[0]])
    d1233 = np.zeros([maskedImage_1d.shape[0]])

    d1311 = np.zeros([maskedImage_1d.shape[0]])
    d1312 = np.zeros([maskedImage_1d.shape[0]])
    d1313 = np.zeros([maskedImage_1d.shape[0]])
    d1321 = np.zeros([maskedImage_1d.shape[0]])
    d1322 = np.zeros([maskedImage_1d.shape[0]])
    d1323 = np.zeros([maskedImage_1d.shape[0]])
    d1331 = np.zeros([maskedImage_1d.shape[0]])
    d1332 = np.zeros([maskedImage_1d.shape[0]])
    d1333 = np.zeros([maskedImage_1d.shape[0]])

    d2111 = np.zeros([maskedImage_1d.shape[0]])
    d2112 = np.zeros([maskedImage_1d.shape[0]])
    d2113 = np.zeros([maskedImage_1d.shape[0]])
    d2121 = np.zeros([maskedImage_1d.shape[0]])
    d2122 = np.zeros([maskedImage_1d.shape[0]])
    d2123 = np.zeros([maskedImage_1d.shape[0]])
    d2131 = np.zeros([maskedImage_1d.shape[0]])
    d2132 = np.zeros([maskedImage_1d.shape[0]])
    d2133 = np.zeros([maskedImage_1d.shape[0]])

    d2211 = np.zeros([maskedImage_1d.shape[0]])
    d2212 = np.zeros([maskedImage_1d.shape[0]])
    d2213 = np.zeros([maskedImage_1d.shape[0]])
    d2221 = np.zeros([maskedImage_1d.shape[0]])
    d2222 = np.zeros([maskedImage_1d.shape[0]])
    d2223 = np.zeros([maskedImage_1d.shape[0]])
    d2231 = np.zeros([maskedImage_1d.shape[0]])
    d2232 = np.zeros([maskedImage_1d.shape[0]])
    d2233 = np.zeros([maskedImage_1d.shape[0]])

    d2311 = np.zeros([maskedImage_1d.shape[0]])
    d2312 = np.zeros([maskedImage_1d.shape[0]])
    d2313 = np.zeros([maskedImage_1d.shape[0]])
    d2321 = np.zeros([maskedImage_1d.shape[0]])
    d2322 = np.zeros([maskedImage_1d.shape[0]])
    d2323 = np.zeros([maskedImage_1d.shape[0]])
    d2331 = np.zeros([maskedImage_1d.shape[0]])
    d2332 = np.zeros([maskedImage_1d.shape[0]])
    d2333 = np.zeros([maskedImage_1d.shape[0]])

    d3111 = np.zeros([maskedImage_1d.shape[0]])
    d3112 = np.zeros([maskedImage_1d.shape[0]])
    d3113 = np.zeros([maskedImage_1d.shape[0]])
    d3121 = np.zeros([maskedImage_1d.shape[0]])
    d3122 = np.zeros([maskedImage_1d.shape[0]])
    d3123 = np.zeros([maskedImage_1d.shape[0]])
    d3131 = np.zeros([maskedImage_1d.shape[0]])
    d3132 = np.zeros([maskedImage_1d.shape[0]])
    d3133 = np.zeros([maskedImage_1d.shape[0]])

    d3211 = np.zeros([maskedImage_1d.shape[0]])
    d3212 = np.zeros([maskedImage_1d.shape[0]])
    d3213 = np.zeros([maskedImage_1d.shape[0]])
    d3221 = np.zeros([maskedImage_1d.shape[0]])
    d3222 = np.zeros([maskedImage_1d.shape[0]])
    d3223 = np.zeros([maskedImage_1d.shape[0]])
    d3231 = np.zeros([maskedImage_1d.shape[0]])
    d3232 = np.zeros([maskedImage_1d.shape[0]])
    d3233 = np.zeros([maskedImage_1d.shape[0]])

    d3311 = np.zeros([maskedImage_1d.shape[0]])
    d3312 = np.zeros([maskedImage_1d.shape[0]])
    d3313 = np.zeros([maskedImage_1d.shape[0]])
    d3321 = np.zeros([maskedImage_1d.shape[0]])
    d3322 = np.zeros([maskedImage_1d.shape[0]])
    d3323 = np.zeros([maskedImage_1d.shape[0]])
    d3331 = np.zeros([maskedImage_1d.shape[0]])
    d3332 = np.zeros([maskedImage_1d.shape[0]])
    d3333 = np.zeros([maskedImage_1d.shape[0]])

    t11 = np.zeros([maskedImage_1d.shape[0]])
    t12 = np.zeros([maskedImage_1d.shape[0]])
    t13 = np.zeros([maskedImage_1d.shape[0]])
    t21 = np.zeros([maskedImage_1d.shape[0]])
    t22 = np.zeros([maskedImage_1d.shape[0]])
    t23 = np.zeros([maskedImage_1d.shape[0]])
    t31 = np.zeros([maskedImage_1d.shape[0]])
    t32 = np.zeros([maskedImage_1d.shape[0]])
    t33 = np.zeros([maskedImage_1d.shape[0]])

    dervXConv = np.zeros([maskedImage_1d.shape[0]])
    dervYConv = np.zeros([maskedImage_1d.shape[0]])
    dervZConv = np.zeros([maskedImage_1d.shape[0]])

    # Second Derivatives

    derXX = np.zeros([maskedImage_1d.shape[0]])
    derYY = np.zeros([maskedImage_1d.shape[0]])
    derZZ = np.zeros([maskedImage_1d.shape[0]])
    derX = np.zeros([maskedImage_1d.shape[0]])
    derXY = np.zeros([maskedImage_1d.shape[0]])
    derXZ = np.zeros([maskedImage_1d.shape[0]])
    derY = np.zeros([maskedImage_1d.shape[0]])
    derYZ = np.zeros([maskedImage_1d.shape[0]])

    tempImgArrayCurr = maskedImage_1d
    tempImgArrayPrev = maskedImage_1d

    gausKernal = ndimage.gaussian_filter1d([1.0, 2.0, 3.0], 1)

    l2_norm = np.linalg.norm(image_1d - maskedImage_1d)
    print("L2 norm = ", l2_norm)

    N = num_of_steps
    # while (l2_norm > tol)
    for n in range(0, N):
        outConvX = convolution3DX(presentImage, masked_width,
                                  masked_height, masked_depth, 1, gausKernal, kernelSize)
        outConvXY = convolution3DY(outConvX, masked_width, masked_height,
                                   masked_depth, 1, gausKernal, kernelSize)
        outConvXYZ = convolution3DZ(outConvXY, masked_width, masked_height,
                                    masked_depth, 1, gausKernal, kernelSize)

        # First derivatives
        dervXConv = derivative3DX(outConvXYZ, masked_width, masked_height,
                                  masked_depth, stepSize)
        dervYConv = derivative3DY(outConvXYZ, masked_width, masked_height,
                                  masked_depth, stepSize)
        dervZConv = derivative3DZ(outConvXYZ, masked_width, masked_height,
                                  masked_depth, stepSize)

        # Second Derivatives

        derXX = derivative3DXX(image, masked_width, masked_height,
                               masked_depth, stepSize)
        derYY = derivative3DYY(image, masked_width, masked_height,
                               masked_depth, stepSize)
        derZZ = derivative3DZZ(image, masked_width, masked_height,
                               masked_depth, stepSize)
        derX = derivative3DX(image, masked_width, masked_height,
                             masked_depth, stepSize)
        derXY = derivative3DY(derX, masked_width, masked_height,
                              masked_depth, stepSize)
        derXZ = derivative3DZ(derX, masked_width, masked_height,
                              masked_depth, stepSize)
        derY = derivative3DY(image, masked_width, masked_height,
                             masked_depth, stepSize)
        derYZ = derivative3DZ(derY, masked_width, masked_height,
                              masked_depth, stepSize)

        for i in range(0, masked_width * masked_height * masked_depth):
            if n == 0:
                tempImgArrayPrev[i] = image[i]

            tempImgArrayCurr[i] = image[i]

            norm_i_square = dervXConv[i] ** 2 + dervYConv[i] ** 2 + dervZConv[i] ** 2
            normi = np.sqrt(norm_i_square)
            xy_norm_i_square = dervXConv[i] ** 2 + dervYConv[i] ** 2
            xy_norm_i = np.sqrt(xy_norm_i_square)

            v1 = np.zeros(3)
            v2 = np.zeros(3)
            v3 = np.zeros(3)
            e1 = np.zeros([3, 3])
            e2 = np.zeros([3, 3])
            e3 = np.zeros([3, 3])
            e4 = np.zeros([3, 3])
            e5 = np.zeros([3, 3])
            e6 = np.zeros([3, 3])

            if norm_i_square == 0:
                v1 = [1, 0, 0]
                v2 = [0, 1, 0]
                v3 = [0, 0, 1]
            else:
                if xy_norm_i_square == 0:
                    v1 = [0, 0, dervZConv[i] / normi]
                    v2 = [0, 0, 0]
                    v3 = [0, 0, 0]
                else:
                    v1 = [dervXConv[i] / normi, dervYConv[i] / normi, dervZConv[i] / normi]
                    v2 = [dervYConv[i] / xy_norm_i, dervXConv[i] / xy_norm_i, 0]
                    v3 = [dervXConv[i] * dervZConv[i] / xy_norm_i * normi,
                          dervYConv[i] * dervZConv[i] / xy_norm_i * normi,
                          -xy_norm_i_square / xy_norm_i * normi]

            m1 = charbonnier_diff(norm_i_square, 0.561268)
            m2 = 1
            m3 = 1
            m4 = np.sqrt(m1 * m2)
            m5 = np.sqrt(m1 * m3)
            m6 = 1

            e1 = [(v1[0] ** 2, v1[0] * v1[1], v1[0] * v1[2]),
                  (v1[1] * v1[0], v1[1] ** 2, v1[1] * v1[2]),
                  (v1[2] * v1[0], v1[2] * v1[1], v1[2] ** 2)]
            e2 = [(v2[0] ** 2, v2[0] * v2[1], v2[0] * v2[2]),
                  (v2[1] * v2[0], v2[1] ** 2, v2[1] * v2[2]),
                  (v2[2] * v2[0], v2[2] * v2[1], v2[2] ** 2)]
            e3 = [(v3[0] ** 2, v3[0] * v3[1], v3[0] * v3[2]),
                  (v3[1] * v3[0], v3[1] ** 2, v3[1] * v3[2]),
                  (v3[2] * v3[0], v3[2] * v3[1], v3[2] ** 2)]

            e4 = [((np.sqrt(2) * v1[0] * v2[0]),
                   (v1[1] * v2[0] + v2[1] * v1[0]) / np.sqrt(2), (v1[2] * v2[0] + v2[2] * v1[0]) / np.sqrt(2)),
                  ((v1[0] * v2[1] + v2[0] * v1[1]) / np.sqrt(2),
                   (np.sqrt(2) * v1[1] * v2[1]), (v1[2] * v2[1] + v2[2] * v1[1]) / np.sqrt(2)),
                  ((v1[0] * v2[2] + v2[0] * v1[2]) / np.sqrt(2),
                   (v1[1] * v2[2] + v2[1] * v1[2]) / np.sqrt(2), (np.sqrt(2) * v1[2] * v2[2]))]

            e5 = [((np.sqrt(2) * v1[0] * v3[0]),
                   (v1[1] * v3[0] + v3[1] * v1[0]) / np.sqrt(2), (v1[2] * v3[0] + v3[2] * v1[0]) / np.sqrt(2)),
                  ((v1[0] * v3[1] + v3[0] * v1[1]) / np.sqrt(2),
                   (np.sqrt(2) * v1[1] * v3[1]), (v1[2] * v3[1] + v3[2] * v1[1]) / np.sqrt(2)),
                  ((v1[0] * v3[2] + v3[0] * v1[2]) / np.sqrt(2), (v1[1] * v3[2] + v3[1] * v1[2]) / np.sqrt(2),
                   (np.sqrt(2) * v1[2] * v3[2]))]

            e6 = [((np.sqrt(2) * v2[0] * v3[0]), (v2[1] * v3[0] + v3[1] * v2[0]) / np.sqrt(2),
                   (v2[2] * v3[0] + v3[2] * v2[0]) / np.sqrt(2)),
                  ((v2[0] * v3[1] + v3[0] * v2[1]) / np.sqrt(2),
                   (np.sqrt(2) * v2[1] * v3[1]), (v2[2] * v3[1] + v3[2] * v2[1]) / np.sqrt(2)),
                  ((v2[0] * v3[2] + v3[0] * v2[2]) / np.sqrt(2),
                   (v2[1] * v3[2] + v3[1] * v2[2]) / np.sqrt(2), (np.sqrt(2) * v2[2] * v3[2]))]

            d1111[i] = m1 * e1[0][0] * e1[0][0] + m2 * e2[0][0] * e2[0][0] + \
                       m3 * e3[0][0] * e3[0][0] + m4 * e4[0][0] * \
                       e4[0][0] + m5 * e5[0][0] * e5[0][0] + m6 * e6[0][0] * e6[0][0]
            d1112[i] = m1 * e1[0][0] * e1[0][1] + m2 * e2[0][0] * e2[0][1] + \
                       m3 * e3[0][0] * e3[0][1] + m4 * e4[0][0] * \
                       e4[0][1] + m5 * e5[0][0] * e5[0][1] + m6 * e6[0][0] * e6[0][1]
            d1113[i] = m1 * e1[0][0] * e1[0][2] + m2 * e2[0][0] * e2[0][2] + \
                       m3 * e3[0][0] * e3[0][2] + m4 * e4[0][0] * \
                       e4[0][2] + m5 * e5[0][0] * e5[0][2] + m6 * e6[0][0] * e6[0][2]
            d1121[i] = m1 * e1[0][0] * e1[1][0] + m2 * e2[0][0] * e2[1][0] + \
                       m3 * e3[0][0] * e3[1][0] + m4 * e4[0][0] * \
                       e4[1][0] + m5 * e5[0][0] * e5[1][0] + m6 * e6[0][0] * e6[1][0]
            d1122[i] = m1 * e1[0][0] * e1[1][1] + m2 * e2[0][0] * e2[1][1] + \
                       m3 * e3[0][0] * e3[1][1] + m4 * e4[0][0] * \
                       e4[1][1] + m5 * e5[0][0] * e5[1][1] + m6 * e6[0][0] * e6[1][1]
            d1123[i] = m1 * e1[0][0] * e1[1][2] + m2 * e2[0][0] * e2[1][2] + \
                       m3 * e3[0][0] * e3[1][2] + m4 * e4[0][0] * \
                       e4[1][2] + m5 * e5[0][0] * e5[1][2] + m6 * e6[0][0] * e6[1][2]
            d1131[i] = m1 * e1[0][0] * e1[2][0] + m2 * e2[0][0] * e2[2][0] + \
                       m3 * e3[0][0] * e3[2][0] + m4 * e4[0][0] * \
                       e4[2][0] + m5 * e5[0][0] * e5[2][0] + m6 * e6[0][0] * e6[2][0]
            d1132[i] = m1 * e1[0][0] * e1[2][1] + m2 * e2[0][0] * e2[2][1] + \
                       m3 * e3[0][0] * e3[2][1] + m4 * e4[0][0] * \
                       e4[2][1] + m5 * e5[0][0] * e5[2][1] + m6 * e6[0][0] * e6[2][1]
            d1133[i] = m1 * e1[0][0] * e1[2][2] + m2 * e2[0][0] * e2[2][2] + \
                       m3 * e3[0][0] * e3[2][2] + m4 * e4[0][0] * \
                       e4[2][2] + m5 * e5[0][0] * e5[2][2] + m6 * e6[0][0] * e6[2][2]

            d1211[i] = m1 * e1[0][1] * e1[0][0] + m2 * e2[0][1] * e2[0][0] + \
                       m3 * e3[0][1] * e3[0][0] + m4 * e4[0][1] * \
                       e4[0][0] + m5 * e5[0][1] * e5[0][0] + m6 * e6[0][1] * e6[0][0]
            d1212[i] = m1 * e1[0][1] * e1[0][1] + m2 * e2[0][1] * e2[0][1] + \
                       m3 * e3[0][1] * e3[0][1] + m4 * e4[0][1] * \
                       e4[0][1] + m5 * e5[0][1] * e5[0][1] + m6 * e6[0][1] * e6[0][1]
            d1213[i] = m1 * e1[0][1] * e1[0][2] + m2 * e2[0][1] * e2[0][2] + \
                       m3 * e3[0][1] * e3[0][2] + m4 * e4[0][1] * \
                       e4[0][2] + m5 * e5[0][1] * e5[0][2] + m6 * e6[0][1] * e6[0][2]
            d1221[i] = m1 * e1[0][1] * e1[1][0] + m2 * e2[0][1] * e2[1][0] + \
                       m3 * e3[0][1] * e3[1][0] + m4 * e4[0][1] * \
                       e4[1][0] + m5 * e5[0][1] * e5[1][0] + m6 * e6[0][1] * e6[1][0]
            d1222[i] = m1 * e1[0][1] * e1[1][1] + m2 * e2[0][1] * e2[1][1] + \
                       m3 * e3[0][1] * e3[1][1] + m4 * e4[0][1] * \
                       e4[1][1] + m5 * e5[0][1] * e5[1][1] + m6 * e6[0][1] * e6[1][1]
            d1223[i] = m1 * e1[0][1] * e1[1][2] + m2 * e2[0][1] * e2[1][2] + \
                       m3 * e3[0][1] * e3[1][2] + m4 * e4[0][1] * \
                       e4[1][2] + m5 * e5[0][1] * e5[1][2] + m6 * e6[0][1] * e6[1][2]
            d1231[i] = m1 * e1[0][1] * e1[2][0] + m2 * e2[0][1] * e2[2][0] + \
                       m3 * e3[0][1] * e3[2][0] + m4 * e4[0][1] * \
                       e4[2][0] + m5 * e5[0][1] * e5[2][0] + m6 * e6[0][1] * e6[2][0]
            d1232[i] = m1 * e1[0][1] * e1[2][1] + m2 * e2[0][1] * e2[2][1] + \
                       m3 * e3[0][1] * e3[2][1] + m4 * e4[0][1] * \
                       e4[2][1] + m5 * e5[0][1] * e5[2][1] + m6 * e6[0][1] * e6[2][1]
            d1233[i] = m1 * e1[0][1] * e1[2][2] + m2 * e2[0][1] * e2[2][2] + \
                       m3 * e3[0][1] * e3[2][2] + m4 * e4[0][1] * \
                       e4[2][2] + m5 * e5[0][1] * e5[2][2] + m6 * e6[0][1] * e6[2][2]

            d1311[i] = m1 * e1[0][2] * e1[0][0] + m2 * e2[0][2] * e2[0][0] + \
                       m3 * e3[0][2] * e3[0][0] + m4 * e4[0][2] * \
                       e4[0][0] + m5 * e5[0][2] * e5[0][0] + m6 * e6[0][2] * e6[0][0]
            d1312[i] = m1 * e1[0][2] * e1[0][1] + m2 * e2[0][2] * e2[0][1] + \
                       m3 * e3[0][2] * e3[0][1] + m4 * e4[0][2] * \
                       e4[0][1] + m5 * e5[0][2] * e5[0][1] + m6 * e6[0][2] * e6[0][1]
            d1313[i] = m1 * e1[0][2] * e1[0][2] + m2 * e2[0][2] * e2[0][2] + \
                       m3 * e3[0][2] * e3[0][2] + m4 * e4[0][2] * \
                       e4[0][2] + m5 * e5[0][2] * e5[0][2] + m6 * e6[0][2] * e6[0][2]
            d1321[i] = m1 * e1[0][2] * e1[1][0] + m2 * e2[0][2] * e2[1][0] + \
                       m3 * e3[0][2] * e3[1][0] + m4 * e4[0][2] * \
                       e4[1][0] + m5 * e5[0][2] * e5[1][0] + m6 * e6[0][2] * e6[1][0]
            d1322[i] = m1 * e1[0][2] * e1[1][1] + m2 * e2[0][2] * e2[1][1] + \
                       m3 * e3[0][2] * e3[1][1] + m4 * e4[0][2] * \
                       e4[1][1] + m5 * e5[0][2] * e5[1][1] + m6 * e6[0][2] * e6[1][1]
            d1323[i] = m1 * e1[0][2] * e1[1][2] + m2 * e2[0][2] * e2[1][2] + \
                       m3 * e3[0][2] * e3[1][2] + m4 * e4[0][2] * \
                       e4[1][2] + m5 * e5[0][2] * e5[1][2] + m6 * e6[0][2] * e6[1][2]
            d1331[i] = m1 * e1[0][2] * e1[2][0] + m2 * e2[0][2] * e2[2][0] + \
                       m3 * e3[0][2] * e3[2][0] + m4 * e4[0][2] * \
                       e4[2][0] + m5 * e5[0][2] * e5[2][0] + m6 * e6[0][2] * e6[2][0]
            d1332[i] = m1 * e1[0][2] * e1[2][1] + m2 * e2[0][2] * e2[2][1] + \
                       m3 * e3[0][2] * e3[2][1] + m4 * e4[0][2] * \
                       e4[2][1] + m5 * e5[0][2] * e5[2][1] + m6 * e6[0][2] * e6[2][1]
            d1333[i] = m1 * e1[0][2] * e1[2][2] + m2 * e2[0][2] * e2[2][2] + \
                       m3 * e3[0][2] * e3[2][2] + m4 * e4[0][2] * \
                       e4[2][2] + m5 * e5[0][2] * e5[2][2] + m6 * e6[0][2] * e6[2][2]

            d2111[i] = m1 * e1[1][0] * e1[0][0] + m2 * e2[1][0] * e2[0][0] + \
                       m3 * e3[1][0] * e3[0][0] + m4 * e4[1][0] * \
                       e4[0][0] + m5 * e5[1][0] * e5[0][0] + m6 * e6[1][0] * e6[0][0]
            d2112[i] = m1 * e1[1][0] * e1[0][1] + m2 * e2[1][0] * e2[0][1] + \
                       m3 * e3[1][0] * e3[0][1] + m4 * e4[1][0] * \
                       e4[0][1] + m5 * e5[1][0] * e5[0][1] + m6 * e6[1][0] * e6[0][1]
            d2113[i] = m1 * e1[1][0] * e1[0][2] + m2 * e2[1][0] * e2[0][2] + \
                       m3 * e3[1][0] * e3[0][2] + m4 * e4[1][0] * \
                       e4[0][2] + m5 * e5[1][0] * e5[0][2] + m6 * e6[1][0] * e6[0][2]
            d2121[i] = m1 * e1[1][0] * e1[1][0] + m2 * e2[1][0] * e2[1][0] + \
                       m3 * e3[1][0] * e3[1][0] + m4 * e4[1][0] * \
                       e4[1][0] + m5 * e5[1][0] * e5[1][0] + m6 * e6[1][0] * e6[1][0]
            d2122[i] = m1 * e1[1][0] * e1[1][1] + m2 * e2[1][0] * e2[1][1] + \
                       m3 * e3[1][0] * e3[1][1] + m4 * e4[1][0] * \
                       e4[1][1] + m5 * e5[1][0] * e5[1][1] + m6 * e6[1][0] * e6[1][1]
            d2123[i] = m1 * e1[1][0] * e1[1][2] + m2 * e2[1][0] * e2[1][2] + \
                       m3 * e3[1][0] * e3[1][2] + m4 * e4[1][0] * \
                       e4[1][2] + m5 * e5[1][0] * e5[1][2] + m6 * e6[1][0] * e6[1][2]
            d2131[i] = m1 * e1[1][0] * e1[2][0] + m2 * e2[1][0] * e2[2][0] + \
                       m3 * e3[1][0] * e3[2][0] + m4 * e4[1][0] * \
                       e4[2][0] + m5 * e5[1][0] * e5[2][0] + m6 * e6[1][0] * e6[2][0]
            d2132[i] = m1 * e1[1][0] * e1[2][1] + m2 * e2[1][0] * e2[2][1] + \
                       m3 * e3[1][0] * e3[2][1] + m4 * e4[1][0] * \
                       e4[2][1] + m5 * e5[1][0] * e5[2][1] + m6 * e6[1][0] * e6[2][1]
            d2133[i] = m1 * e1[1][0] * e1[2][2] + m2 * e2[1][0] * e2[2][2] + \
                       m3 * e3[1][0] * e3[2][2] + m4 * e4[1][0] * \
                       e4[2][2] + m5 * e5[1][0] * e5[2][2] + m6 * e6[1][0] * e6[2][2]

            d2211[i] = m1 * e1[1][1] * e1[0][0] + m2 * e2[1][1] * e2[0][0] + \
                       m3 * e3[1][1] * e3[0][0] + m4 * e4[1][1] * \
                       e4[0][0] + m5 * e5[1][1] * e5[0][0] + m6 * e6[1][1] * e6[0][0]
            d2212[i] = m1 * e1[1][1] * e1[0][1] + m2 * e2[1][1] * e2[0][1] + \
                       m3 * e3[1][1] * e3[0][1] + m4 * e4[1][1] * \
                       e4[0][1] + m5 * e5[1][1] * e5[0][1] + m6 * e6[1][1] * e6[0][1]
            d2213[i] = m1 * e1[1][1] * e1[0][2] + m2 * e2[1][1] * e2[0][2] + \
                       m3 * e3[1][1] * e3[0][2] + m4 * e4[1][1] * \
                       e4[0][2] + m5 * e5[1][1] * e5[0][2] + m6 * e6[1][1] * e6[0][2]
            d2221[i] = m1 * e1[1][1] * e1[1][0] + m2 * e2[1][1] * e2[1][0] + \
                       m3 * e3[1][1] * e3[1][0] + m4 * e4[1][1] * \
                       e4[1][0] + m5 * e5[1][1] * e5[1][0] + m6 * e6[1][1] * e6[1][0]
            d2222[i] = m1 * e1[1][1] * e1[1][1] + m2 * e2[1][1] * e2[1][1] + \
                       m3 * e3[1][1] * e3[1][1] + m4 * e4[1][1] * \
                       e4[1][1] + m5 * e5[1][1] * e5[1][1] + m6 * e6[1][1] * e6[1][1]
            d2223[i] = m1 * e1[1][1] * e1[1][2] + m2 * e2[1][1] * e2[1][2] + \
                       m3 * e3[1][1] * e3[1][2] + m4 * e4[1][1] * \
                       e4[1][2] + m5 * e5[1][1] * e5[1][2] + m6 * e6[1][1] * e6[1][2]
            d2231[i] = m1 * e1[1][1] * e1[2][0] + m2 * e2[1][1] * e2[2][0] + \
                       m3 * e3[1][1] * e3[2][0] + m4 * e4[1][1] * \
                       e4[2][0] + m5 * e5[1][1] * e5[2][0] + m6 * e6[1][1] * e6[2][0]
            d2232[i] = m1 * e1[1][1] * e1[2][1] + m2 * e2[1][1] * e2[2][1] + \
                       m3 * e3[1][1] * e3[2][1] + m4 * e4[1][1] * \
                       e4[2][1] + m5 * e5[1][1] * e5[2][1] + m6 * e6[1][1] * e6[2][1]
            d2233[i] = m1 * e1[1][1] * e1[2][2] + m2 * e2[1][1] * e2[2][2] + \
                       m3 * e3[1][1] * e3[2][2] + m4 * e4[1][1] * \
                       e4[2][2] + m5 * e5[1][1] * e5[2][2] + m6 * e6[1][1] * e6[2][2]

            d2311[i] = m1 * e1[1][2] * e1[0][0] + m2 * e2[1][2] * e2[0][0] + \
                       m3 * e3[1][2] * e3[0][0] + m4 * e4[1][2] * \
                       e4[0][0] + m5 * e5[1][2] * e5[0][0] + m6 * e6[1][2] * e6[0][0]
            d2312[i] = m1 * e1[1][2] * e1[0][1] + m2 * e2[1][2] * e2[0][1] + \
                       m3 * e3[1][2] * e3[0][1] + m4 * e4[1][2] * \
                       e4[0][1] + m5 * e5[1][2] * e5[0][1] + m6 * e6[1][2] * e6[0][1]
            d2313[i] = m1 * e1[1][2] * e1[0][2] + m2 * e2[1][2] * e2[0][2] + \
                       m3 * e3[1][2] * e3[0][2] + m4 * e4[1][2] * \
                       e4[0][2] + m5 * e5[1][2] * e5[0][2] + m6 * e6[1][2] * e6[0][2]
            d2321[i] = m1 * e1[1][2] * e1[1][0] + m2 * e2[1][2] * e2[1][0] + \
                       m3 * e3[1][2] * e3[1][0] + m4 * e4[1][2] * \
                       e4[1][0] + m5 * e5[1][2] * e5[1][0] + m6 * e6[1][2] * e6[1][0]
            d2322[i] = m1 * e1[1][2] * e1[1][1] + m2 * e2[1][2] * e2[1][1] + \
                       m3 * e3[1][2] * e3[1][1] + m4 * e4[1][2] * \
                       e4[1][1] + m5 * e5[1][2] * e5[1][1] + m6 * e6[1][2] * e6[1][1]
            d2323[i] = m1 * e1[1][2] * e1[1][2] + m2 * e2[1][2] * e2[1][2] + \
                       m3 * e3[1][2] * e3[1][2] + m4 * e4[1][2] * \
                       e4[1][2] + m5 * e5[1][2] * e5[1][2] + m6 * e6[1][2] * e6[1][2]
            d2331[i] = m1 * e1[1][2] * e1[2][0] + m2 * e2[1][2] * e2[2][0] + \
                       m3 * e3[1][2] * e3[2][0] + m4 * e4[1][2] * \
                       e4[2][0] + m5 * e5[1][2] * e5[2][0] + m6 * e6[1][2] * e6[2][0]
            d2332[i] = m1 * e1[1][2] * e1[2][1] + m2 * e2[1][2] * e2[2][1] + \
                       m3 * e3[1][2] * e3[2][1] + m4 * e4[1][2] * \
                       e4[2][1] + m5 * e5[1][2] * e5[2][1] + m6 * e6[1][2] * e6[2][1]
            d2333[i] = m1 * e1[1][2] * e1[2][2] + m2 * e2[1][2] * e2[2][2] + \
                       m3 * e3[1][2] * e3[2][2] + m4 * e4[1][2] * \
                       e4[2][2] + m5 * e5[1][2] * e5[2][2] + m6 * e6[1][2] * e6[2][2]

            d3111[i] = m1 * e1[2][0] * e1[0][0] + m2 * e2[2][0] * e2[0][0] + \
                       m3 * e3[2][0] * e3[0][0] + m4 * e4[2][0] * \
                       e4[0][0] + m5 * e5[2][0] * e5[0][0] + m6 * e6[2][0] * e6[0][0]
            d3112[i] = m1 * e1[2][0] * e1[0][1] + m2 * e2[2][0] * e2[0][1] + \
                       m3 * e3[2][0] * e3[0][1] + m4 * e4[2][0] * \
                       e4[0][1] + m5 * e5[2][0] * e5[0][1] + m6 * e6[2][0] * e6[0][1]
            d3113[i] = m1 * e1[2][0] * e1[0][2] + m2 * e2[2][0] * e2[0][2] + \
                       m3 * e3[2][0] * e3[0][2] + m4 * e4[2][0] * \
                       e4[0][2] + m5 * e5[2][0] * e5[0][2] + m6 * e6[2][0] * e6[0][2]
            d3121[i] = m1 * e1[2][0] * e1[1][0] + m2 * e2[2][0] * e2[1][0] + \
                       m3 * e3[2][0] * e3[1][0] + m4 * e4[2][0] * \
                       e4[1][0] + m5 * e5[2][0] * e5[1][0] + m6 * e6[2][0] * e6[1][0]
            d3122[i] = m1 * e1[2][0] * e1[1][1] + m2 * e2[2][0] * e2[1][1] + \
                       m3 * e3[2][0] * e3[1][1] + m4 * e4[2][0] * \
                       e4[1][1] + m5 * e5[2][0] * e5[1][1] + m6 * e6[2][0] * e6[1][1]
            d3123[i] = m1 * e1[2][0] * e1[1][2] + m2 * e2[2][0] * e2[1][2] + \
                       m3 * e3[2][0] * e3[1][2] + m4 * e4[2][0] * \
                       e4[1][2] + m5 * e5[2][0] * e5[1][2] + m6 * e6[2][0] * e6[1][2]
            d3131[i] = m1 * e1[2][0] * e1[2][0] + m2 * e2[2][0] * e2[2][0] + \
                       m3 * e3[2][0] * e3[2][0] + m4 * e4[2][0] * \
                       e4[2][0] + m5 * e5[2][0] * e5[2][0] + m6 * e6[2][0] * e6[2][0]
            d3132[i] = m1 * e1[2][0] * e1[2][1] + m2 * e2[2][0] * e2[2][1] + \
                       m3 * e3[2][0] * e3[2][1] + m4 * e4[2][0] * \
                       e4[2][1] + m5 * e5[2][0] * e5[2][1] + m6 * e6[2][0] * e6[2][1]
            d3133[i] = m1 * e1[2][0] * e1[2][2] + m2 * e2[2][0] * e2[2][2] + \
                       m3 * e3[2][0] * e3[2][2] + m4 * e4[2][0] * \
                       e4[2][2] + m5 * e5[2][0] * e5[2][2] + m6 * e6[2][0] * e6[2][2]

            d3211[i] = m1 * e1[2][1] * e1[0][0] + m2 * e2[2][1] * e2[0][0] + \
                       m3 * e3[2][1] * e3[0][0] + m4 * e4[2][1] * \
                       e4[0][0] + m5 * e5[2][1] * e5[0][0] + m6 * e6[2][1] * e6[0][0]
            d3212[i] = m1 * e1[2][1] * e1[0][1] + m2 * e2[2][1] * e2[0][1] + \
                       m3 * e3[2][1] * e3[0][1] + m4 * e4[2][1] * \
                       e4[0][1] + m5 * e5[2][1] * e5[0][1] + m6 * e6[2][1] * e6[0][1]
            d3213[i] = m1 * e1[2][1] * e1[0][2] + m2 * e2[2][1] * e2[0][2] + \
                       m3 * e3[2][1] * e3[0][2] + m4 * e4[2][1] * \
                       e4[0][2] + m5 * e5[2][1] * e5[0][2] + m6 * e6[2][1] * e6[0][2]
            d3221[i] = m1 * e1[2][1] * e1[1][0] + m2 * e2[2][1] * e2[1][0] + \
                       m3 * e3[2][1] * e3[1][0] + m4 * e4[2][1] * \
                       e4[1][0] + m5 * e5[2][1] * e5[1][0] + m6 * e6[2][1] * e6[1][0]
            d3222[i] = m1 * e1[2][1] * e1[1][1] + m2 * e2[2][1] * e2[1][1] + \
                       m3 * e3[2][1] * e3[1][1] + m4 * e4[2][1] * \
                       e4[1][1] + m5 * e5[2][1] * e5[1][1] + m6 * e6[2][1] * e6[1][1]
            d3223[i] = m1 * e1[2][1] * e1[1][2] + m2 * e2[2][1] * e2[1][2] + \
                       m3 * e3[2][1] * e3[1][2] + m4 * e4[2][1] * \
                       e4[1][2] + m5 * e5[2][1] * e5[1][2] + m6 * e6[2][1] * e6[1][2]
            d3231[i] = m1 * e1[2][1] * e1[2][0] + m2 * e2[2][1] * e2[2][0] + \
                       m3 * e3[2][1] * e3[2][0] + m4 * e4[2][1] * \
                       e4[2][0] + m5 * e5[2][1] * e5[2][0] + m6 * e6[2][1] * e6[2][0]
            d3232[i] = m1 * e1[2][1] * e1[2][1] + m2 * e2[2][1] * e2[2][1] + \
                       m3 * e3[2][1] * e3[2][1] + m4 * e4[2][1] * \
                       e4[2][1] + m5 * e5[2][1] * e5[2][1] + m6 * e6[2][1] * e6[2][1]
            d3233[i] = m1 * e1[2][1] * e1[2][2] + m2 * e2[2][1] * e2[2][2] + \
                       m3 * e3[2][1] * e3[2][2] + m4 * e4[2][1] * \
                       e4[2][2] + m5 * e5[2][1] * e5[2][2] + m6 * e6[2][1] * e6[2][2]

            d3311[i] = m1 * e1[2][2] * e1[0][0] + m2 * e2[2][2] * e2[0][0] + \
                       m3 * e3[2][2] * e3[0][0] + m4 * e4[2][2] * \
                       e4[0][0] + m5 * e5[2][2] * e5[0][0] + m6 * e6[2][2] * e6[0][0]
            d3312[i] = m1 * e1[2][2] * e1[0][1] + m2 * e2[2][2] * e2[0][1] + \
                       m3 * e3[2][2] * e3[0][1] + m4 * e4[2][2] * \
                       e4[0][1] + m5 * e5[2][2] * e5[0][1] + m6 * e6[2][2] * e6[0][1]
            d3313[i] = m1 * e1[2][2] * e1[0][2] + m2 * e2[2][2] * e2[0][2] + \
                       m3 * e3[2][2] * e3[0][2] + m4 * e4[2][2] * \
                       e4[0][2] + m5 * e5[2][2] * e5[0][2] + m6 * e6[2][2] * e6[0][2]
            d3321[i] = m1 * e1[2][2] * e1[1][0] + m2 * e2[2][2] * e2[1][0] + \
                       m3 * e3[2][2] * e3[1][0] + m4 * e4[2][2] * \
                       e4[1][0] + m5 * e5[2][2] * e5[1][0] + m6 * e6[2][2] * e6[1][0]
            d3322[i] = m1 * e1[2][2] * e1[1][1] + m2 * e2[2][2] * e2[1][1] + \
                       m3 * e3[2][2] * e3[1][1] + m4 * e4[2][2] * \
                       e4[1][1] + m5 * e5[2][2] * e5[1][1] + m6 * e6[2][2] * e6[1][1]
            d3323[i] = m1 * e1[2][2] * e1[1][2] + m2 * e2[2][2] * e2[1][2] + \
                       m3 * e3[2][2] * e3[1][2] + m4 * e4[2][2] * \
                       e4[1][2] + m5 * e5[2][2] * e5[1][2] + m6 * e6[2][2] * e6[1][2]
            d3331[i] = m1 * e1[2][2] * e1[2][0] + m2 * e2[2][2] * e2[2][0] + \
                       m3 * e3[2][2] * e3[2][0] + m4 * e4[2][2] * \
                       e4[2][0] + m5 * e5[2][2] * e5[2][0] + m6 * e6[2][2] * e6[2][0]
            d3332[i] = m1 * e1[2][2] * e1[2][1] + m2 * e2[2][2] * e2[2][1] + \
                       m3 * e3[2][2] * e3[2][1] + m4 * e4[2][2] * \
                       e4[2][1] + m5 * e5[2][2] * e5[2][1] + m6 * e6[2][2] * e6[2][1]
            d3333[i] = m1 * e1[2][2] * e1[2][2] + m2 * e2[2][2] * e2[2][2] + \
                       m3 * e3[2][2] * e3[2][2] + m4 * e4[2][2] * \
                       e4[2][2] + m5 * e5[2][2] * e5[2][2] + m6 * e6[2][2] * e6[2][2]

            t11[i] = d1111[i] * derXX[i] + d1112[i] * derXY[i] + \
                     d1113[i] * derXZ[i] + d1121[i] * derXY[i] + d1122[i] * \
                     derYY[i] + d1123[i] * derYZ[i] + d1131[i] * derXZ[i] + \
                     d1132[i] * derYZ[i] + d1133[i] * derZZ[i]
            t12[i] = d1211[i] * derXX[i] + d1212[i] * derXY[i] + d1213[i] * \
                     derXZ[i] + d1221[i] * derXY[i] + d1222[i] * \
                     derYY[i] + d1223[i] * derYZ[i] + d1231[i] * \
                     derXZ[i] + d1232[i] * derYZ[i] + d1233[i] * derZZ[i];
            t13[i] = d1311[i] * derXX[i] + d1312[i] * derXY[i] + \
                     d1313[i] * derXZ[i] + d1321[i] * derXY[i] + d1322[i] * \
                     derYY[i] + d1323[i] * derYZ[i] + d1331[i] * derXZ[i] + \
                     d1332[i] * derYZ[i] + d1333[i] * derZZ[i]

            t21[i] = d2111[i] * derXX[i] + d2112[i] * derXY[i] + \
                     d2113[i] * derXZ[i] + d2121[i] * derXY[i] + d2122[i] * \
                     derYY[i] + d2123[i] * derYZ[i] + d2131[i] * derXZ[i] + \
                     d2132[i] * derYZ[i] + d2133[i] * derZZ[i]
            t22[i] = d2211[i] * derXX[i] + d2212[i] * derXY[i] + \
                     d2213[i] * derXZ[i] + d2221[i] * derXY[i] + d2222[i] * \
                     derYY[i] + d2223[i] * derYZ[i] + d2231[i] * derXZ[i] + \
                     d2232[i] * derYZ[i] + d2233[i] * derZZ[i]
            t23[i] = d2311[i] * derXX[i] + d2312[i] * derXY[i] + d2313[i] * \
                     derXZ[i] + d2321[i] * derXY[i] + d2322[i] * \
                     derYY[i] + d2323[i] * derYZ[i] + d2331[i] * \
                     derXZ[i] + d2332[i] * derYZ[i] + d2333[i] * derZZ[i]

            t31[i] = d3111[i] * derXX[i] + d3112[i] * derXY[i] + \
                     d3113[i] * derXZ[i] + d3121[i] * derXY[i] + d3122[i] * \
                     derYY[i] + d3123[i] * derYZ[i] + d3131[i] * derXZ[i] + \
                     d3132[i] * derYZ[i] + d3133[i] * derZZ[i]
            t32[i] = d3211[i] * derXX[i] + d3212[i] * derXY[i] + \
                     d3213[i] * derXZ[i] + d3221[i] * derXY[i] + d3222[i] * \
                     derYY[i] + d3223[i] * derYZ[i] + d3231[i] * derXZ[i] + \
                     d3232[i] * derYZ[i] + d3233[i] * derZZ[i]
            t33[i] = d3311[i] * derXX[i] + d3312[i] * derXY[i] + d3313[i] * \
                     derXZ[i] + d3321[i] * derXY[i] + d3322[i] * \
                     derYY[i] + d3323[i] * derYZ[i] + d3331[i] * \
                     derXZ[i] + d3332[i] * derYZ[i] + d3333[i] * derZZ[i]

    dervXXD11 = derivative3DXX(t11, masked_width, masked_height, masked_depth, stepSize)
    dervYYD22 = derivative3DYY(t22, masked_width, masked_height, masked_depth, stepSize)
    dervZZD33 = derivative3DZZ(t33, masked_width, masked_height, masked_depth, stepSize)
    dervYD21 = derivative3DX(t21, masked_width, masked_height, masked_depth, stepSize)
    dervXYD21 = derivative3DY(dervYD21, masked_width, masked_height, masked_depth, stepSize)
    dervZD31 = derivative3DX(t31, masked_width, masked_height, masked_depth, stepSize)
    dervXZD31 = derivative3DZ(dervZD31, masked_width, masked_height, masked_depth, stepSize)
    dervXD12 = derivative3DY(t12, masked_width, masked_height, masked_depth, stepSize)
    dervYXD12 = derivative3DX(dervXD12, masked_width, masked_height, masked_depth, stepSize)
    dervZD32 = derivative3DY(t32, masked_width, masked_height, masked_depth, stepSize)
    dervYZD32 = derivative3DZ(dervZD32, masked_width, masked_height, masked_depth, stepSize)
    dervXD13 = derivative3DZ(t13, masked_width, masked_height, masked_depth, stepSize)
    dervZXD13 = derivative3DX(dervXD13, masked_width, masked_height, masked_depth, stepSize)
    dervYD23 = derivative3DZ(t23, masked_width, masked_height, masked_depth, stepSize)
    dervZYD23 = derivative3DY(dervYD23, masked_width, masked_height, masked_depth, stepSize)

    alpha = (4 * n + 2) / (2 * n + 3)

    randArrTraceIndex = 0

    for i in range (0, masked_width*masked_height*masked_depth):
        if i == randomPixels[randArrTraceIndex]:
            randArrTraceIndex+=1

        presentImage[i] = alpha * (presentImage[i] - timeStep * (
                    dervXXD11[i] + dervXYD21[i] + dervXZD31[i] + dervYXD12[i] +
                    dervYYD22[i] + dervYZD32[i] + dervZXD13[i] + dervZYD23[i] +
                    dervZZD33[i])) + (1 - alpha) * tempImgArrayPrev[i]

        tempImgArrayPrev[i] = tempImgArrayCurr[i]

    print("Run Time = ", datetime.now() - currentTime)
    display_image(presentImage)

    return 0
