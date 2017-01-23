from setuptools import setup, find_packages
from distutils.core import setup, Extension
import numpy.distutils.misc_util
import sys, os

here = os.path.abspath(os.path.dirname(__file__))
README = open(os.path.join(here, 'README.txt')).read()
NEWS = open(os.path.join(here, 'NEWS.txt')).read()

version = '1.0'

# http://packages.python.org/distribute/setuptools.html#declaring-dependencies
install_requires = [
    # 
]

include_dirs = numpy.distutils.misc_util.get_numpy_include_dirs()
include_dirs.append('../../include')

setup(name='albopictus',
    version=version,
    description="Environmentally-driven population dynamics model of Aedes albopictus",
    long_description=README + '\n\n' + NEWS,
    # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.'
    ],
    keywords=['stage','age','structured','gridded','global','diapause','breeding','egg','larva','pupa','adult','mosquito','temperature','precipitation','density','photoperiod','survival','development','fecundity','Bayesian','difference','daily'],
    author='Kamil Erguler',
    author_email='k.erguler@cyi.ac.cy',
    url = 'https://github.com/kerguler/albopictus',
    download_url = "https://github.com/kerguler/albopictus/tarball/%s" %(version),
    license='GPLv3',
    packages=find_packages('src'),
    package_dir = {'': 'src'},
    include_package_data=True,
    package_data={'albopictus': ['data/*.json']},
    zip_safe=False,
    install_requires=install_requires,
    py_modules=['albopictus/__init__'],
    ext_modules=[
        Extension("albopictus.modelAalbopictus03", ["src/albopictus/incubator03.c", "src/albopictus/modelAalbopictus03.c"]),
        Extension("albopictus.modelAalbopictus", ["src/albopictus/gamma.c", "src/albopictus/incubator.c", "src/albopictus/modelAalbopictus.c"]),
        Extension("albopictus.modelStochCHIKV", ["src/albopictus/ran_gen.c", "src/albopictus/spop.c", "src/albopictus/gamma.c", "src/albopictus/modelStochCHIKV.c"]),
        Extension("albopictus.modelStochSand", ["src/albopictus/ran_gen.c", "src/albopictus/spop.c", "src/albopictus/gamma.c", "src/albopictus/modelStochSand.c"])],
    include_dirs=include_dirs
)
