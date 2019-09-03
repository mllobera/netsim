

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
These notebooks provide examples of how to use the different functionality offered in **NetSim**

Feature
-------

**NetSim** uses a *influence weighted distance transform* or *iwdt* a purposefully designed
algorithm that allows for the use of an '*influence layer*'. This layer can be used to generate
the convergence of posterior paths onto earlier paths hence producing path networks that tend to be
more simple and 'natural' looking.

Installation
============

The simplest way to install this package is by cloning the git repository and following these *recommended*
steps:

1. Download and install `miniconda <https://conda.io/projects/conda/en/latest/user-guide/install/index.html?highlight=conda>`_

2. Generate a new environment using the ``netsim.yml``

   >>> conda env create -f netsim. yml

3. Activate the new environment,

   >>> activate netsim

4. Use ``pip`` to install the **NetSim** package into the ``site-packages`` folder 
   in the *netsim* environment you newly created.

   >>> pip install dist/netsim*.tar.gz


Citation
========

This work is currently under review.  For the time being please cite this work as follows:

    Llobera, M. 2019. Netsim package (Version 0.0.1). https://github.com/mllobera/netsim

  








