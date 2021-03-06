import setuptools
from os.path import join, sep, split, abspath
from glob import glob

with open("README.md", "r") as fh:
    long_description = fh.read()

here = split(abspath(__file__))[0]
datapath = here + sep + 'share' + sep + 'lofarcabling' + sep + 'layouts'
datafiles = glob(datapath + sep + '*')

setuptools.setup(
    name="lofarcabling",
    version="1.0",
    author="Thomas Tijsma",
    author_email="thomas.tijsma@gmail.com",
    description="A project that designs cable layouts for LOFAR",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lofar-astron/lofar-cabling",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    zip_safe = False,
    data_files = [('share/lofarcabling/layouts',
                  datafiles)]
)
