# tinyHPIPM 
<!-- a stripped down version of the HPIPM solver -->
This project is an opinionated and stripped down version of [HPIPM](https://github.com/giaf/hpipm) packaged with the linear algebra subroutines in [BLASFEO](https://github.com/giaf/blasfeo).
I created it based on my own personal needs (which were mostly oriented around QCQPs and the Python interface), but if you want me to add something, don't hesitate to open a PR.  
Currently, the main differences with the original project(s) are:
- for BLASFEO: a lot of the compilation options have been removed to keep only those actually used in HPIPM (high performance implementation with panel major matrix format, no external BLAS interface).
- for HPIPM: no tree OCP interface
- the compilation targets for both have been merged, simplified and renamed (based on CPU capabilities instead of vendor architectures)
- no more Make build system, only modern CMake that simplifies the usage from other project
- no more Matlab interface 
- some small bugfixes in the OCP QCQP interface

In the future I would like to add:
- a more extensive Python interface 
- more testing and code analysis to ensure space grade software quality
    - unit testing (probably with CTest in the beginning but I could switch to a more complete C framework)
    - high MCDC coverage
    - use static analysis tools like SonarLint, IKOS, COBRA to ensure standards
- more documentation

