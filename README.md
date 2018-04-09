# Symbolic and Analysis Control package (SAC)
SAC is a set of routines for analysis and control of nonlinear time-delay systems (altough it also works for systems without delays).  It is based on  [MAXIMA](http://maxima.sourceforge.net), a GPL computer algebra system, available on all major platforms.

# Installation
To use this software, you will need a recent version of MAXIMA.  SAC has been developed using version 5.37.2 and GCL, but it should work with any version supporting structure data types.

To install it, just copy all the files from [https://github.com/alxmz2/sac.git](https://github.com/alxmz2/sac.git) to a directory, and either 
* start maxima from that directory, or
* set the enviroment variable MAXIMA_USERDIR to include the directory containing the software.  On linux, this can be done by adding the line 

> export MAXIMA_USERDIR=$MAXIMA_USERDIR:/path/to/sac

to the .bashrc file located in your home directory.

Also, you can unpack the documentation file in a directory (e.g. /path/to/sac/doc )
and open the index.html file with your browser.  Some features require to have javascript enabled.  This documentation includes several examples, as well as a General Hierarchy Diagram, which shows all the defined functions and their dependencies.

# Usage

After starting maxima, just load the main file:

> (%i1) load("sac.mc")$

All the routines will be ready for use. See the documentation for examples.

# About
Most of this version has been rewritten from scratch. The previous one can be found [here](https://sourceforge.net/projects/sac/).

# Credits
The original project SAC was written by

* Luis Alejandro MARQUEZ-MARTINEZ
* Araceli GARATE-GARCIA
* José Ricardo CUESTA-GARCIA
* Eduardo GARCIA-RAMIREZ


