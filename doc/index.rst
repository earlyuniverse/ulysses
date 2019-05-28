.. ulysses documentation master file, created by
   sphinx-quickstart on Fri May 17 16:40:52 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ulysses's documentation!
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Documentation for the Code
**************************


ULSBase
=========

This is the base class.

.. autoclass:: ulysses.ULSBase
    :members:

    .. automethod:: __init__
    .. automethod:: __call__


ulysses.tools
=============
.. automodule:: ulysses.tools
    :members:

Built-in models
===============

.. autoclass:: ulysses.EtaB_1DME
    :members:

    .. automethod:: RHS
    .. automethod:: EtaB

.. autoclass:: ulysses.EtaB_2DME
    :members:

    .. automethod:: RHS
    .. automethod:: EtaB

.. autoclass:: ulysses.EtaB_3DME
    :members:

    .. automethod:: RHS
    .. automethod:: EtaB

.. autoclass:: ulysses.EtaB_1BE
    :members:

    .. automethod:: RHS
    .. automethod:: EtaB

.. autoclass:: ulysses.EtaB_2BE
    :members:

    .. automethod:: RHS
    .. automethod:: EtaB

.. autoclass:: ulysses.EtaB_2Resonant
    :members:

    .. automethod:: RHS
    .. automethod:: EtaB

.. autoclass:: ulysses.EtaB_3DME_Scattering
    :members:

    .. automethod:: RHS
    .. automethod:: EtaB

.. autoclass:: ulysses.EtaB_3DS_Scattering_RHtaur
    :members:

    .. automethod:: RHS
    .. automethod:: EtaB

