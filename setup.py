import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

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
)
