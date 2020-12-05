# OLSApy

***

## OLSApy: Orthogonal Linear Separation Analysis in Python

***
* OLSA is an analysis method of omics data to decompose the complex effects of a perturbagen into basic components.  
* OLSApy is a package for OLSA in python.  
* OLSA can be applied to any kinds of omics data such as RNA-seq, proteome, and so on.  

***

## Dependency
* python >= 3.6
* requirements: numpy, pandas, scipy

***

## Setup
```pip install olsapy```  
For conda, please use conda skeleton  

***

## Usage

    
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
    

***

## Licence
This software is released under the MIT License, see LICENSE.  

***

## Authors
Setsuo Kinoshita, Shotaro Maedera, and Tadahaya Mizuno  

***

## References
* Sci Rep. 2019 Feb 12;9(1):1824. doi: 10.1038/s41598-019-38528-4.  
* http://www.ilincs.org/ilincs/  

***

## Bug Report
If you would like to report any bugs about olsapy, don't hesitate to create an issue on github here, or email me:  
* tadahaya@gmail.com  
* tadahaya@mol.f.u-tokyo.ac.jp  