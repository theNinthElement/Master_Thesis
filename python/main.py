# #from linear_homogeneous_inpainting import linear_homogeneous_FSI
#
# #import cython
# #from EED_3D_inpainting import EED_FSI
#
# from example_add import test
# #from imgCodecs import newExample as m
#
# def print_hi(name):
#     # Use a breakpoint in the code line below to debug your script.
#     print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.
#
#
# # Press the green button in the gutter to run the script.
# if __name__ == '__main__':
#     # linear_homogeneous_FSI()
#     m.test()


import rect
# In zoo.cpp we expose hello() function, and it now exists in the zoo module.
x0, y0, x1, y1 = 1, 2, 3, 4
# zoo.hello is a callable.
rect_obj = rect.PyRectangle()
# Call the C++ hello() function from Python.
rect_obj._cinit_(x0, y0, x1, y1)
print(rect_obj.get_size())
print(rect_obj.get_area())