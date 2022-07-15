# import os
# import pathlib
# from distutils.core import setup
# import pkg_resources
# import setuptools

# VERSION = "0.1"


# def build_packages(base_dir, name_base):
#     """
#     recusively find all the folders and treat them as packages
#     """
#     arr = [name_base]
#     for fname in os.listdir(base_dir):
#         if os.path.isdir(base_dir + fname):
#             """
#             ignore the hidden files
#             """
#             if fname[0] == ".":
#                 continue
#             name = "{}.{}".format(name_base, fname)
#             recursion = build_packages(base_dir + fname + "/", name)
#             if len(recursion) != 0:
#                 [arr.append(rec) for rec in recursion]
#     return arr


# packages = build_packages("sosat-optics/", "sosat-optics")

# setup(
#     name="sosat_optics",
#     version=VERSION,
#     description="Optical simulation of the Simons Observatory Small Aperture Telescope.",
#     author="Grace E. Chesmore, UChicago Lab",
#     author_email="chesmore@uchicago.edu",
#     packages = ['sosat_optics'],
#     package_dir = {'sosat_optics':'sosat-optics'},
# )


# Install requirements
# with pathlib.Path("requirements.txt").open() as requirements_txt:
#     install_requires = [
#         str(requirement)
#         for requirement in pkg_resources.parse_requirements(requirements_txt)
#     ]

# setuptools.setup(
#     install_requires=install_requires,
# )


import os
from setuptools import setup, Command, find_packages

def readme():
    with open('README.md') as f:
        return f.read()

class CleanCommand(Command):
    """Custom clean command to tidy up the project root."""
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        os.system('rm -vrf build dist *.pyc *.tgz *.egg-info')


setup(name = 'sosat_optics',
      version = '0.1',
      description = '''Code for Modeling Polarization in Optics''',
      long_description = readme(),
      author = 'Katie Harrington',
      author_email = 'katie.megan.harrington@gmail.com',
      license = 'MIT',
      packages = ['sosat_optics'],
      package_dir = {'sosat_optics':'sosat-optics'},
      cmdclass={'clean':CleanCommand,},
      #scripts = []
      )