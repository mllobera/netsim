

.. |NetSim_Logo| image:: ../images/Netsim_Logo.png 
   :width: 50px                                       
   :height: 50px                    
   :alt: NetSim                     


|NetSim_Logo|  **NetSim** is a python package for the simulation of the evolution of a path network.

Description
===========

The **NetSim** package was designed with landscape archaeology in mind where it is very common
to group sites into simple groups. These groups might reflect coarse temporal, functional or other
differences (associated with one or combination of features of the sites). These groupings often 
translate into some form of hierarchy that will affect the evolution of a path network.

This package is accompanied by a series of notebooks that can be found in the ``/notebooks`` folder.
These notebooks provide examples of how to use the different functionality offered in **NetSim**.

Feature
-------

**NetSim** uses a *influence weighted distance transform* or *iwdt* a purposefully designed
algorithm that allows for the use of an '*influence layer*'. This layer can be used to generate
the convergence of posterior paths onto earlier paths hence producing path networks that tend to be
more simple and 'natural' looking.

Installation
============

This package has yet to be uploaded to PyPI. There are several ways to install the **NetSim** package.

To run this package you will need to have several python packages installed. The best way to
achieve this is by generating a 'virtual' environment where all the required packages are installed.
This can be easily achieved by `installing miniconda`_ first.

Currently, certain parts of the package need to be compiled on installation. To achieve this, 
you might need to download a compiler into your machine if one is not currently installed. The
following instructions will guide you to achieve this if you are installing the package on a Win10 
machine. *Pending* (Mac and Linux instructions).

For windows
-----------

Installing Microsoft compiler on you machine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the latest on compilers you may want to check `here <https://wiki.python.org/moin/WindowsCompilers>`_ 

You need to install a standalone version of *Visual C++ 14.2* compiler. n.b. You do not need to install
*Visual Studio 2019*.

- Go to `Microsoft Build Tools for Visual Studio 2019 <https://www.visualstudio.com/downloads/#build-tools-for-visual-studio-2019>`_
- Look for the section 'Build Tools for Visual Studio' and install *C++ build tools* . You can use the default
  but make sure that you have the latest versions of:

  - MSVCv142 - VS 2019 C++ x64/x86 build tools
  - Windows 10 SDK

Installing miniconda
^^^^^^^^^^^^^^^^^^^^

Download and install `miniconda <https://conda.io/projects/conda/en/latest/user-guide/install/index.html?highlight=conda>`_
for Python 3.7 or above.

- Create a project folder

- Download the ``netsim.yml`` file into this folder by clicking on the file in the git repo. Click on **raw** button 
  and right-click to `save as...` the file onto your project folder (make sure that you save it with only the '.yml'
  extension)

- Generate a new environment using the ``netsim.yml`` file

  >>> conda env create -f netsim.yml

- Activate the new environment,

   >>> activate netsim

Install NetSim package
^^^^^^^^^^^^^^^^^^^^^^

- Navigate and download the **NetSim** package source file into your project folder

- Use ``pip`` to install the **NetSim** package into the ``site-packages`` folder 
  in the *netsim* environment you newly created.

   >>> pip install dist/netsim*.tar.gz

Documentation
=============

To learn how to use **NetSim** consult the following `documentation <https://netsim.readthedocs.io/>`_


Citation
========

This work is currently under review.  For the time being please cite this work as follows:

    Llobera, M. 2019. Netsim package (Version 0.0.1). https://github.com/mllobera/netsim


License
=======

**NetSim** is an open source `BSD-license <../../../LICENSE.rst>`_ package 










