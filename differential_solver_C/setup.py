from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np
from pathlib import Path

this_dir = Path(__file__).parent
readme = (this_dir / "README.md").read_text(encoding="utf-8")

ext = Extension(
    name="SIR_diffusion_wrap",
    sources=["SIR_diffusion_wrap.pyx", "SIR_diffusion.c"],
    include_dirs=[np.get_include(), "."],
    extra_compile_args=["-O3", "-fPIC"],
    language="c",
)

setup(
    name="sir-diffusion",                     # distribution name on PyPI
    version="0.1.0",                          # bump this when you change code
    description="Cython-accelerated SIR diffusion solver",
    long_description=readme,
    long_description_content_type="text/markdown",
    author="Your Name",
    author_email="you@example.com",
    url="https://example.com/your-project-or-github",
    ext_modules=cythonize([ext], compiler_directives={"language_level": 3}),
    python_requires=">=3.8",
    install_requires=["numpy"],              # runtime deps
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: C",
        "Programming Language :: Python :: Implementation :: CPython",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)

