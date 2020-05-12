from setuptools import setup, find_packages


setup(
    name="rebin1d",
    version="0.0.1",
    author="Jae-Joon Lee",
    author_email="lee.j.joon@gmail.com",
    description="A simple 1d rebin package",
    package_dir={'rebin1d': 'rebin1d'},
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
