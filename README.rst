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
* python 3.6
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

 from olsapy import olsa as ol
   
3. generate a DataClass object as follows:

::

 dat = ol.DataClass()

4. load the prepared data file into the generated object as follows:

::

 dat.load(<a path for the data file>)

5. run OLSA and obtain a Result object as follows:

::

 res = ol.olsa(dat)

6. export each result as csv files as follows:

::

 res.export()

7. each result can be extracted as a dataframe if necessary as follows:

::

 dataframe = res.rsm()

* a sample code for running OLSA described below:

::

 from olsapy import olsa as ol
   
 filein = '<file path>'

 #run OLSA simply
 dat = ol.DataClass() #generate a DataClass object
 dat.load(filein) #load data
 res = ol.olsa(dat) #run OLSA and obtain a Result object
 res.export() #save data

 #run OLSA with some options
 df = res.rsm() #.rsm(), etc. extract stored data in a Result object as a dataframe
 dat2 = ol.DataClass()
 dat2.load_df(df) #load dataframe into a DataClass object
 res2 = ol.olsa(dat2,accumulation=0.5) #accumulation determines the vectors subjected to varimax rotation
 res2.export(CM=True,TS=False) #results to be exported can be chosen.

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