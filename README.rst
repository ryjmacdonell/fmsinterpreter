FMSinterpreter
==============

Analysis and plotting scripts for FMS90 and FMSpy output files.

Installation
------------
To add to your local Python packages, clone the repository and use setup.py
to install via::

    $ git clone https://github.com/ryjmacdonell/fmsinterpreter.git
    $ python setup.py install

This will also install several command line scripts located in ``bin/``. In
particular, the ``getdefault`` script followed a script name will generate
an input file with the default inputs.

Usage
-----
The command line scripts will readily interpret FMS data based on common
tasks such as plotting adiabatic populations or showing trends in spawn
geometries. For more advanced tasks, FMSinterpreter modules can be
incorporated into Python scripts via::

    from fmsinterpreter import [module]

where ``[module]`` is the module for the required task.

Requirements
------------
Requires the `GeomTools library <https://github.com/ryjmacdonell/geomtools>`_
and at least Python 3.3, NumPy v1.6.0, SciPy v0.9.0
