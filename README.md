# 🤏 tinyHPIPM
<!-- a stripped down version of the HPIPM solver -->
This project is an opinionated and stripped down version of [`HPIPM`](https://github.com/giaf/hpipm) packaged with the linear algebra subroutines in [`BLASFEO`](https://github.com/giaf/blasfeo).
I created it based on my own personal needs, that were mostly oriented around QCQPs, Python, CMake and ARM processors at the time I started this project.

## Main differences with the old projects

- *no more compilation targets*: they were a pain to figure out, to configure, to maintain and to use. What's more, on ARM platforms the compiler actually did enough auto-vectorizaion to match the performances using the generic computational kernels written in C. So I decided to ditch the ~800k lines of assembly and only keep the C kernels.
- *no more Matlab interface*: I never used it and it was a pain to maintain. I might add it back in the future if I need it.
- *no more tree OCP interface*: I never used it and it was a pain to maintain. I might add it back in the future if I need it.
- *only CMake build system and no more Make*: to be able to easily integrate `tinyHPIPM` in other projects and to add more comprehensive tooling and unit testing, I decided to only keep a modern CMake build system based on good practices.
- *no more single precision versions*: first of all, they weren't useful for my applications. Second of all, the old genericity system (see [this issue](https://github.com/tudoroancea/tinyhpipm/issues/2)) was very hard to maintain, mainly due to excessive use of macros and consequent poor `clangd` analysis and IDE integration. I therefore chose to only keep the double precision versions.
- *modern C code style*: I tried to make the code more readable and more in line with modern C code style.
- *only the BLAS/LAPACK operations that we need*: `BLASFEO` was very comprehensive but a lot of code was never used in `HPIPM` and refactoring it is very time consuming

## Roadmap

### Features
- [ ] extensive unit testing, code coverage and static analysis (using tools like [SonarLint](https://www.sonarsource.com/products/sonarlint/), [IKOS](https://github.com/NASA-SW-VnV/ikos), [COBRA](https://github.com/nimble-code/Cobra)) to ensure space grade software quality.
- [ ] extensive documentation and examples including:
    - the data storage format in `mat`
    - the BLAS/LAPACK subroutines implemented
    - the solver interfaces (dense and OCP)
- [ ] simplify the interface even further and remove the QP interfaces to only keep the QCQP ones (since they comprises the former with little to no overhead).
- [ ] revamp the Python interface to make it more _pythonic_ (already started in [this file](python/tinyhpipm/new_ocp_qcqp.py)).

### Refactors
- [ ] change all integer and floating point data types to generic types `tinyhpipm_int_t` and `tinyhpipm_float_t` defined in `tinyhpipm/common.h` to allow for easy porting to other types
