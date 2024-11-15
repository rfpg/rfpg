

### Makefile.pdlibbuilder ###

Helper makefile for Pure Data external libraries.
Written by Katja Vetter March-June 2015 for the public domain. No warranties.
Inspired by Hans Christoph Steiner's Makefile Template and Stephan Beal's
ShakeNMake.

GNU make version >= 3.81 required.


### characteristics ###


* defines build settings based on autodetected OS and architecture
* defines rules to build Pd class- or lib executables from C or C++ sources
* defines rules for libdir installation
* defines convenience targets for developer and user
* evaluates implicit dependencies for non-clean builds


### basic usage ###


In your Makefile, define your Pd lib name and class files, and include
Makefile.pdlibbuilder at the end of the Makefile. Like so:


      # Makefile for mylib

      lib.name = mylib

      class.sources = myclass1.c myclass2.c

      datafiles = myclass1-help.pd myclass2-help.pd README.txt LICENSE.txt

      include Makefile.pdlibbuilder


For files in class.sources it is assumed that class basename == source file
basename. The default target builds all classes as individual executables
with Pd's default extension for the platform. For anything more than the
most basic usage, read the documentation sections in Makefile.pdlibbuilder.


### examples ###


Here is one deployment example of the Makefile.pdlibbuilder approach:

http://sourceforge.net/p/pure-data/svn/HEAD/tree/trunk/externals/miXed/cyclone/Makefile.cyclone

More examples will be referenced here when they are available. 
