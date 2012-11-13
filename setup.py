# setuptools installation of Hollow
# setup.py written by Oliver Beckstein [2010]
# ... and cast into the Public Domain...

from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, find_packages

# Dynamically calculate the version based on VERSION.
version = __import__('hollow').__version__


setup(name="Hollow",
      version=version,
      description="Volume filling of protein structures",
      long_description="""\ 
Hollow helps to visualize and quantify
voids, pockets, and channels in protein structures.

It generates fake atoms that fill space inside a structure specified
in the PDB format. The output is a standard PDB file containing the
fake atoms. The molecular surface of these atoms is an excellent
approximation to the cavity and can be easily rendered in standard
molecular visualization software such as PyMOL, VMD or Chimera.

In complements programs such as HOLE, CAVER, MOLE, SPHGEN, and CASTp.
""",
      author="Franz Gruswitz and Bosco Ho",
      author_email="boscoh@gmail.com",
      license="unknown",  # probably GPLv2
      url="http://hollow.sourceforge.net/",
      keywords="science molecular structure",
      packages=find_packages(exclude=['scripts']),
      package_data = {'hollow': ['*.txt', 
                                 'pdbstruct/*.txt', 'pdbstruct/template.pdb']},
      scripts = ["scripts/hollow",
                 "scripts/hollow-asa",
                 "scripts/hollow-volume",
                 ],
      #install_requires=['psyco'],  # psyco is really optional...
      zip_safe = True,
)

      
