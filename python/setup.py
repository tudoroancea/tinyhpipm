from setuptools import setup, find_packages

setup(
    name="tinyhpipm",
    version="0.0.0",
    description="Python interface to tinyHPIPM",
    url="http://github.com/tudoroancea/tinyhpipm",
    author="Tudor Oancea - Andrea Zanelli - Gianluca Frison",
    license="BSD-2",
    packages=find_packages(),
    zip_safe=False,
    install_requires=["numpy"],
)
