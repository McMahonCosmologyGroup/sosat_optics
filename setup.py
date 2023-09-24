import os
from setuptools import setup, find_packages

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
    version="0.1.1",
    description="Optical simulation of the Simons Observatory Small Aperture Telescope.",
    long_description=open("README.md", "r", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    author="Grace E. Chesmore, Alex Thomas, UChicago Lab",
    author_email="chesmore@uchicago.edu, agthomas@uchicago.edu",
    packages=packages,
    package_dir={"sosat_optics": "sosat_optics"},
    install_requires=["numpy", "matplotlib", "scipy==1.10.1"]
)