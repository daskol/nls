#!/usr/bin/env python2

"""NLS: Non-Linear Schrodinger equation solver.

NLS is scientific package that provides ability to solve effeffectively
non-linear schrodinger equation with reservoir. These equation describes
exciton-polariton condensation in microcavities. NLS is built on native
fortran code and is based on certain natural abstraction layer that wraps
native solver. These features are reason that makes calcualtions with NLS
fast.
"""

from numpy.distutils.core import setup
from numpy.distutils.core import Extension as FortranExtension

DOCLINES = (__doc__ or '').split('\n')

CLASSIFIERS = """\
Development Status :: 4 - Beta
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: MIT License
Programming Language :: Fortran
Programming Language :: Python
Programming Language :: Python :: 2
Programming Language :: Python :: 2.6
Programming Language :: Python :: 2.7
Programming Language :: Python :: Implementation :: CPython
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

PLATFORMS = [
    'Windows',
    'Linux',
    'Solaris',
    'Mac OS-X',
    'Unix'
]

MAJOR = 0
MINOR = 1
PATCH = 0

VERSION = '{0:d}.{1:d}.{2:d}'.format(MAJOR, MINOR, PATCH)


def setup_package():
    setup(name='nls',
         version=VERSION,
         description = DOCLINES[0],
         long_description = '\n'.join(DOCLINES[2:]),
         url='https://github.com/daskol/nls',
         author='Daniel Bershatsky',
         author_email='daniel.bershatsky@skolkovotech.ru',
         license='MIT',
         platforms=PLATFORMS,
         classifiers=[line for line in CLASSIFIERS.split('\n') if line],
         packages=[
            'nls'
         ],
         ext_modules=[
            FortranExtension(
                name='nls.native',
                sources=[
                    'nls/nls.f95'
                ],
                libraries=[
                    'blas',
                ],
            ),
         ],
         install_requires=[
            'numpy',
            'scipy',
            'matplotlib',
         ],
         include_package_data=True,
         zip_safe=False
    )


if __name__ == '__main__':
    setup_package()
