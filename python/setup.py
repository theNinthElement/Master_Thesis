from setuptools import setup, Extension
from Cython.Build import cythonize

#os.environ["CC"] = "g++"

setup(ext_modules = cythonize(Extension(
      "c_inpainting",                   # our Cython source
      sources = ["c_inpainting.pyx"],  # additional source file(s)
      language = "c++"             # generate C++ code
  )))