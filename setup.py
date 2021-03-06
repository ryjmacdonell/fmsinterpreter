"""
Setup script for the fmsinterpreter package
"""
from setuptools import setup
from setuptools import find_packages


def readme():
    with open('README.rst') as f:
        return f.read()


setup(
    name='fmsinterpreter',
    version='0.1',
    description=('Tools for analysis of FMS90 and nomad output files'),
    long_description=readme(),
    keywords='fmsinterpreter fms dynamics chemistry',
    url='https://github.com/ryjmacdonell/fmsinterpreter',
    author='Ryan J. MacDonell',
    author_email='rmacd054@uottawa.ca',
    license='MIT',
    packages=find_packages(),
    scripts=[
        'bin/contourplot',
        'bin/denplot',
        'bin/geninput',
        'bin/getdefault',
        'bin/getgeoms',
        'bin/histplot',
        'bin/pesplot',
        'bin/popassign',
        'bin/popplot',
        'bin/scatterplot',
        'bin/trpesplot'
             ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Topic :: Scientific/Engineering :: Chemistry'
                 ],
    install_requires=['numpy>=1.6.0', 'scipy>=0.9.0']
      )
