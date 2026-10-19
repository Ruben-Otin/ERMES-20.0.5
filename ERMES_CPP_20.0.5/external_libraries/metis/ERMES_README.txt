METIS and GKlib bundled with ERMES 20.0.5
==========================================

   metis/   METIS 5 ( https://github.com/KarypisLab/METIS )  - libmetis sources and metis.h
   GKlib/   GKlib   ( https://github.com/KarypisLab/GKlib )  - sources and headers

Both are copyright the Regents of the University of Minnesota, Apache License 2.0
( metis/LICENSE, GKlib/LICENSE.txt; a few GKlib files carry other licenses, listed in
GKlib/LICENSES.md ).

They are compiled automatically with ERMES ( Linux Makefile and Visual Studio project ) and
provide the nested dissection ordering of the internal direct solver. Nothing has to be installed.

Only change with respect to the original sources: metis/include/metis.h starts with

   #define IDXTYPEWIDTH 32
   #define REALTYPEWIDTH 32

( METIS's own build system normally prepends these two lines ).

To update METIS: replace metis/libmetis, metis/include/metis.h ( adding the two lines ) and the
GKlib folders with the new sources from GitHub, and rebuild ERMES.
