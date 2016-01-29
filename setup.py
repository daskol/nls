from setuptools import setup

setup(name='nls',
     version='0.0.0',
     description='Non-linear Schrodinger equation solver.',
     url='https://github.com/daskol/nls',
     author='Daniel Bershatsky',
     author_email='daniel.bershatsky@skolkovotech.ru',
     license='MIT',
     packages=['nls'],
     install_requires=[
        'numpy',
        'matplotlib',
     ],
     include_package_data=True,
     zip_safe=False)