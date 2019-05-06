#/usr/bin/env python

"""
setup.py for snapshot module
"""

from distutils.core import setup, Extension
import sys
if sys.version_info >= (3, 0):
    sys.exit("Python3 not supported.")

snap_base = 'gadget_binary_snap/'
snap_module = Extension(
    '_binarysnap',
    sources=[snap_base + f for f in (
        'snap_wrap.cpp',
        'snapclass.cpp'
    )],
    define_macros=[('SWIG', None),
                   ('OUTPUTPOTENTIALS', None),
                   ('OUTPUTTIMESTEP', None)
    ],
)
    
setup(
    name='pygadgetsnap',
    version='0.1',
    description='''Gadget snapshot (binary, hdf5) wrapper for python.''',
    url='',
    author='Kyle Oman',
    author_email='koman@astro.rug.nl',
    license='',
    ext_modules=[snap_module],
    py_modules=['gadget_hdf5_snap', 'gadget_binary_snap'],
    install_requires=['numpy', 'h5py'],
    include_package_data=True,
    zip_safe=False
)

