FOAM-FSI_oldSchool-thijs
========

Is a branch of original start of Thijs Gillebaart his contributions to FOAM-FSI.
This branch holds all of the code at the start, including unfinished not ideal parts.

See master branch for latest update on this.

First, compile foam-extend-3.1, the nextRelease branch.

To compile the FSI library:

``` bash
cd src/
./Allmake
cd ../applications/solvers
./Allwmake
cd ../applications/utilities
```

Prerequisites
-----------

gcc 4.8 or higher due to C++11 features.
