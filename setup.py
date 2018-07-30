from setuptools import setup, find_packages
from distutils.core import setup, Extension
import numpy.distutils.misc_util
import sys, os

# --------------------------------------------------------------------------
# https://stackoverflow.com/questions/38523941/change-cythons-naming-rules-for-so-files
# Thanks to hoefling for positing a solution to the Cython's naming problem

from distutils.command.install_lib import install_lib as _install_lib

def batch_rename(src, dst, src_dir_fd=None, dst_dir_fd=None):
    '''Same as os.rename, but returns the renaming result.'''
    os.rename(src, dst,
              src_dir_fd=src_dir_fd,
              dst_dir_fd=dst_dir_fd)
    return dst

class _CommandInstallCythonized(_install_lib):
    def __init__(self, *args, **kwargs):
        _install_lib.__init__(self, *args, **kwargs)

    def install(self):
        import re
        # let the distutils' install_lib do the hard work
        outfiles = _install_lib.install(self)
        # batch rename the outfiles:
        # for each file, match string between
        # second last and last dot and trim it
        matcher = re.compile('\.([^.]+)\.so$')
        return [batch_rename(file, re.sub(matcher, '.so', file))
                for file in outfiles]

# --------------------------------------------------------------------------
    
here = os.path.abspath(os.path.dirname(__file__))
README = open(os.path.join(here, 'README.txt')).read()
NEWS = open(os.path.join(here, 'NEWS.txt')).read()

version = '1.9'

# http://packages.python.org/distribute/setuptools.html#declaring-dependencies

include_dirs = numpy.distutils.misc_util.get_numpy_include_dirs()
include_dirs.append('../../include')

setup(name='albopictus',
    version=version,
    description="Large-scale environment-driven population dynamics and disease spread models for vector-borne diseases",
    long_description=README + '\n\n' + NEWS,
    # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.'
    ],
    keywords=['stage','age','structured','gridded','global','diapause','breeding','egg','larva','pupa','adult','mosquito','temperature','precipitation','density','photoperiod','survival','development','fecundity','Bayesian','difference','daily','albopictus','sandfly','chikv','chikungunya','phlebotomus','papatasi'],
    author='Kamil Erguler',
    author_email='k.erguler@cyi.ac.cy',
    url = 'https://github.com/kerguler/albopictus',
    download_url = "https://github.com/kerguler/albopictus/tarball/%s" %(version),
    license='GPLv3',
    cmdclass={
        'install_lib': _CommandInstallCythonized
        },
    packages=find_packages('src'),
    package_dir = {'': 'src'},
    include_package_data=True,
    package_data={'albopictus': ['data/*.json']},
    zip_safe=False,
    install_requires=[],
    py_modules=[
        'albopictus/__init__',
        'albopictus/readModel/__init__',
        'albopictus/setPrior/__init__',
        'albopictus/plotPos/__init__',
        'albopictus/accessory/__init__',
        'albopictus/population/__init__'
        ],
    ext_modules=[
        Extension("albopictus.modelAalbopictus03", ["src/albopictus/incubator03.c", "src/albopictus/modelAalbopictus03.c"]),
        Extension("albopictus.modelAalbopictus08", ["src/albopictus/gamma.c", "src/albopictus/incubator.c", "src/albopictus/modelAalbopictus08.c"]),
        Extension("albopictus.modelAalbopictus13", ["src/albopictus/gamma.c", "src/albopictus/incubator.c", "src/albopictus/modelAalbopictus13.c"]),
        # Extension("albopictus.modelAalbopictus18", ["src/albopictus/gamma.c", "src/albopictus/spop.c", "src/albopictus/modelAalbopictus18.c"]),
        Extension("albopictus.modelStochCHIKV", ["src/albopictus/ran_gen.c", "src/albopictus/spop01.c", "src/albopictus/gamma.c", "src/albopictus/modelStochCHIKV.c"]),
        Extension("albopictus.modelStochSand", ["src/albopictus/ran_gen.c", "src/albopictus/spop01.c", "src/albopictus/gamma.c", "src/albopictus/modelStochSand.c"])
        ],
    include_dirs=include_dirs
)
