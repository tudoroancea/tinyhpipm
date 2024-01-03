# 🤏 tinyHPIPM 
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

> WARNING: This project is in its early stages, a lot can change unexpectedly (especially in the C interface that I have not yet finished cleaning up), so this is currently not production ready.

In the future I would like to add:
- automatic compilation target selection
- a more extensive Python interface 
- more testing and code analysis to ensure space grade software quality
    - unit testing (probably with CTest in the beginning but I could switch to a more complete C framework)
    - high MCDC coverage
    - use static analysis tools like [SonarLint](https://www.sonarsource.com/products/sonarlint/), [IKOS](https://github.com/NASA-SW-VnV/ikos), [COBRA](https://github.com/nimble-code/Cobra) to ensure standards like MISRA-C.
- more documentation


## New compilation targets

The table below is a summary of the correspondance between the old compilation targets of BLASFEO and the new compilation targets of tinyHPIPM.

| old target | new target | cflags | asm flags |
| --- | --- | --- | --- |
| TARGET_X86_AMD_BARCELONA | TARGET_X86_SSE3 | -m32 -msse3 |  |
| TARGET_X86_AMD_JAGUAR | TARGET_X86_AVX | -m32 -mavx |  |
| TARGET_X64_INTEL_CORE | TARGET_X64_SSE3 | -m64 -msse3 |  |
| TARGET_X64_INTEL_SANDY_BRIDGE | TARGET_X64_AVX | -m64 -mavx |  |
| TARGET_X64_AMD_BULLDOZER | TARGET_X64_AVX_FMA | -m64 -mavx -mfma |  |
| TARGET_X64_INTEL_HASWELL | TARGET_X64_AVX2_FMA | -m64 -mavx2 -mfma |  |
| TARGET_X64_INTEL_SKYLAKE_X | TARGET_X64_AVX512_FMA | -m64 -mavx512f -mavx512vl -mfma |  |
| TARGET_ARMV7A_... | TARGET_ARMV7A_NEON_VPFV3 | -marm -mfloat-abi=hard -mfpu=neon-vfpv3 | -mfpu=neon-vfpv4 |
| TARGET_ARMV8A_... | TARGET_ARMV8A_NEON2_VPFV4 | -march=armv8-a+crc+crypto+simd |  |
| TARGET_GENERIC | TARGET_GENERIC |  |  |
