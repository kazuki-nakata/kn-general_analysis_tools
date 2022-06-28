from setuptools import setup, find_packages  # , Extension
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy as np

if __name__ == "__main__":
    extensions = [Extension("*", ["*.pyx"])]

    setup(
        name="knool",
        version="0.1",
        description="This product summarize program codes used in my study",
        author="K. Nakata, N. Tamaru, K. Kai, S. Nakano and D. K. Seguchi",
        author_email="knakata0924@gmail.com",
        url="https://github.com/kazuki-nakata/",
        packages=find_packages(),
        ext_modules=cythonize(["knool/helpers/*.pyx"]),
        include_dirs=[np.get_include()]
        #        install_requires=open('requirements.txt').read().splitlines()
    )
