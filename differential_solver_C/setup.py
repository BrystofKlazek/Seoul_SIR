from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

ext = Extension(
    name="SIR_diffusion_wrap",
    sources=["SIR_diffusion_wrap.pyx", "SIR_diffusion.c"],  # wrapper + your C
    include_dirs=[np.get_include(), "."],
    extra_compile_args=["-O3", "-fPIC"],
    language="c",
)

setup(
    ext_modules=cythonize([ext], compiler_directives={"language_level": 3}),
)

