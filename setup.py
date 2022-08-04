import os
import pathlib
from distutils.core import setup

import pkg_resources
import setuptools

VERSION = "0.1"


def build_packages(base_dir, name_base):
    """
    recusively find all the folders and treat them as packages
    """
    arr = [name_base]
    for fname in os.listdir(base_dir):
        if os.path.isdir(base_dir + fname):
            """
            ignore the hidden files
            """
            if fname[0] == ".":
                continue
            name = "{}.{}".format(name_base, fname)
            recursion = build_packages(base_dir + fname + "/", name)
            if len(recursion) != 0:
                [arr.append(rec) for rec in recursion]
    return arr


packages = build_packages("sosat_optics/", "sosat_optics")

setup(
    name="sosat_optics",
    version=VERSION,
    description="Optical simulation of the Simons Observatory Small Aperture Telescope.",
    author="Grace E. Chesmore, UChicago Lab",
    author_email="chesmore@uchicago.edu",
    packages=["sosat_optics"],
    package_dir={"sosat_optics": "sosat_optics"},
)
