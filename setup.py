import os
from setuptools import setup, find_namespace_packages

setup(
    name="sosat_optics",
    version="0.1.1",
    description="Optical simulation of the Simons Observatory Small Aperture Telescope.",
    long_description=open("README.md", "r", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    author="Grace E. Chesmore, Alex Thomas, UChicago Lab",
    author_email="chesmore@uchicago.edu, agthomas@uchicago.edu",
    packages=["sosat_optics"],
    package_dir={"sosat_optics": "sosat_optics"},
    install_requires=open("requirements.txt").read().splitlines()
)