from setuptools import find_packages

from numpy.distutils.core import setup, Extension

ext1 = Extension(name='lib_fixed_turbation.soilcprocessesf90.turbation',
                 sources=['src/lib_fixed_turbation/soilCProcesses.f90'],
                )



setup(name="lib_fixed_turbation",
      version="0.0.1",
      package_dir={"": "src"},
      packages=find_packages(where="src"),
      ext_modules=[ext1])
