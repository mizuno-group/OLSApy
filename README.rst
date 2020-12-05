========
OLSAPY
========

OLSAPY: Orthogonal Linear Separation Analysis in Python
=======================================================
* OLSA is an analysis method of omics data to decompose the complex effects of a perturbagen into basic components.
* OLSAPY is a package for OLSA in python.
* OLSA can be applied to any kinds of omics data such as RNA-seq, proteome, and so on.

Dependency
=======================================================
* python >= 3.6
* requirements: numpy, pandas, scipy

Setup
=======================================================
::

 pip install olsapy

Usage
=======================================================
1. prepare a profile matrix with variables in rows and samples in columns as a csv file
2. import necessary modules as follows:

::

 from olsapy import olsa
   
3. generate an OLSA object as follows:

::

 dat = olsa.OLSA()

4. load the prepared dataframe (feature x sample matrix) into the generated object as follows:

::

 dat.load_df(<prepared dataframe>)

5. run OLSA as follows:

::

 dat.olsa(dat)

6. get the results as follows:

::

 rsm,rvm,ts,contribution = dat.get_res()

7. visualize the results as follows:

::
 dat.plot(focus=[<sample names of interest>])

* a sample code for running OLSA described below:

::

 from olsapy import olsa
   
 df = <a dataframe of interest (feature x sample)>

 ### calculation
 dat = olsa.OLSA()
 dat.load_df(df)
 dat.olsa(accumulation=0.6) # accumulation determines the vectors subjected to varimax rotation
 rsm,rvm,ts,contribution = dat.get_res()
  
 ### visualization
 dat.plot(focus=[<sample names of interest>])
 
 ### visualization of already prepared data
 dat.plot(rsm,focus=[<sample names of interest>])


Licence
=======================================================
This software is released under the MIT License, see LICENSE.

Authors
=======================================================
Setsuo Kinoshita, Shotaro Maedera, and Tadahaya Mizuno

References
=======================================================
http://www.ilincs.org/ilincs/

Bug Report
=======================================================
If you would like to report any bugs about olsapy, don't hesitate to create an issue on github here, or email me: tadahaya@gmail.com