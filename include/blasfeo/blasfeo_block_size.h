// NOTE: this file is not actually useful for this stripped down version,
// it was mostly for the Column Major format

#ifndef BLASFEO_BLOCK_SIZE_H_
#define BLASFEO_BLOCK_SIZE_H_

#define D_EL_SIZE 8  // double precision
#define S_EL_SIZE 4  // single precision

#if defined(TARGET_X64_AVX512_FMA)
// common
#define CACHE_LINE_SIZE 64  // data cache size: 64 bytes
// double
#define D_PS 8  // panel size
#define D_PLD 8  // 4 // GCD of panel length
#define D_KC 128  // 256 // 192
#define D_NC 144  // 72 //96 //72 // 120 // 512
#define D_MC 2400  // 6000
// single
#define S_PS 16  // panel size
#define S_PLD 4  // GCD of panel length TODO probably 16 when writing assebly
#define S_KC 128  // 256
#define S_NC 128  // 144
#define S_MC 3000

#elif defined(TARGET_X64_AVX2_FMA)
// common
#define CACHE_LINE_SIZE 64  // data cache size: 64 bytes
// double
#define D_PS 4  // panel size
#define D_PLD 4  // 2 // GCD of panel length
#define D_KC 256  // 192
#define D_NC 64  // 96 //72 // 120 // 512
#define D_MC 1500
// single
#define S_PS 8  // panel size
#define S_PLD 4  // 2 // GCD of panel length
#define S_KC 256
#define S_NC 144
#define S_MC 3000

#elif defined(TARGET_X64_AVX)
// common
#define CACHE_LINE_SIZE 64  // data cache size: 64 bytes
// double
#define D_PS 4  // panel size
#define D_PLD 4  // 2 // GCD of panel length
#define D_KC 256  // 320 //256 //320
#define D_NC 72  // 64 //72 //60 // 120
#define D_MC 1000  // 800
// single
#define S_PS 8  // panel size
#define S_PLD 4  // 2 // GCD of panel length
#define S_KC 256
#define S_NC 144
#define S_MC 2000

#elif defined(TARGET_X64_INTEL_CORE)
// common
#define CACHE_LINE_SIZE 64
// double
#define D_PS 4  // panel size
#define D_PLD 4  // 2 // GCD of panel length
#define D_KC 256
#define D_NC 128  // TODO these are just dummy
#define D_MC 3000  // TODO these are just dummy
// single
#define S_PS 4
#define S_PLD 4  // 2
#define S_KC 256
#define S_NC 128  // TODO these are just dummy
#define S_MC 3000  // TODO these are just dummy

#elif defined(TARGET_X64_AVX_FMA)
// common
#define CACHE_LINE_SIZE 64
// double
#define D_PS 4  // panel size
#define D_PLD 4  // 2 // GCD of panel length
#define D_KC 256
#define D_NC 128  // TODO these are just dummy
#define D_MC 3000  // TODO these are just dummy
// single
#define S_PS 4
#define S_PLD 4  // 2
#define S_KC 256
#define S_NC 128  // TODO these are just dummy
#define S_MC 3000  // TODO these are just dummy


#elif defined(TARGET_X86_AVX)
// common
#define CACHE_LINE_SIZE 64
// double
#define D_PS 4  // panel size
#define D_PLD 4  // 2 // GCD of panel length
#define D_KC 256
#define D_NC 128  // TODO these are just dummy
#define D_MC 3000  // TODO these are just dummy
// single
#define S_PS 4
#define S_PLD 4  // 2
#define S_KC 256
#define S_NC 128  // TODO these are just dummy
#define S_MC 3000  // TODO these are just dummy


#elif defined(TARGET_X86_SSE3)
// common
#define CACHE_LINE_SIZE 64
// double
#define D_PS 4  // panel size
#define D_PLD 4  // 2 // GCD of panel length
#define D_KC 256
#define D_NC 128  // TODO these are just dummy
#define D_MC 3000  // TODO these are just dummy
// single
#define S_PS 4
#define S_PLD 4  // 2
#define S_KC 256
#define S_NC 128  // TODO these are just dummy
#define S_MC 3000  // TODO these are just dummy

#elif defined(TARGET_ARMV8A_NEON2_VPFV4)  // should at least work with Apple silicon and ARM cortex A76, for others check below the old code
// common
#define CACHE_LINE_SIZE 64
// double
#define D_PS 4  // panel size
#define D_PLD 4  // 2 // GCD of panel length
#define D_KC 512  // 256
#define D_NC 128  // 256
#define D_MC 6000
// single
#define S_PS 4
#define S_PLD 4  // 2
#define S_KC 512
#define S_NC 256
#define S_MC 6000
#elif defined(TARGET_ARMV7A_NEON_VFPV3)
// common
#define CACHE_LINE_SIZE 64
// double
#define D_PS 4  // panel size
#define D_PLD 4  // 2 // GCD of panel length
#define D_KC 256
#define D_NC 128  // TODO these are just dummy
#define D_MC 3000  // TODO these are just dummy
// single
#define S_PS 4
#define S_PLD 4  // 2
#define S_KC 256
#define S_NC 128  // TODO these are just dummy
#define S_MC 3000  // TODO these are just dummy
#elif defined(TARGET_GENERIC)
// common
#define CACHE_LINE_SIZE 64
#define L1_CACHE_SIZE (32 * 1024)  // L1 data cache size: 32 kB
// double
#define D_PS 4  // panel size
#define D_PLD 4  // 2 // GCD of panel length
#define D_M_KERNEL 4  // max kernel size
#define D_N_KERNEL 4  // max kernel size
#define D_KC 256
#define D_NC 128  // TODO these are just dummy
#define D_MC 3000  // TODO these are just dummy

// single
#define S_PS 4
#define S_PLD 4  // 2
#define S_M_KERNEL 4  // max kernel size
#define S_N_KERNEL 4  // max kernel size
#define S_KC 256
#define S_NC 128  // TODO these are just dummy
#define S_MC 3000  // TODO these are just dummy

// old ARM targets

// #elif defined(TARGET_ARMV8A_APPLE_M1)
// // common
// #define CACHE_LINE_SIZE 64
// // double
// #define D_PS 4  // panel size
// #define D_PLD 4  // 2 // GCD of panel length
// #define D_KC 512  // 256
// #define D_NC 128  // 256
// #define D_MC 6000
// // single
// #define S_PS 4
// #define S_PLD 4  // 2
// #define S_KC 512
// #define S_NC 256
// #define S_MC 6000
// #elif defined(TARGET_ARMV8A_ARM_CORTEX_A76)
// // common
// #define CACHE_LINE_SIZE 64
// // double
// #define D_PS 4  // panel size
// #define D_PLD 4  // 2 // GCD of panel length
// #define D_KC 512  // 256
// #define D_NC 128  // 256
// #define D_MC 6000
// // single
// #define S_PS 4
// #define S_PLD 4  // 2
// #define S_KC 512
// #define S_NC 256
// #define S_MC 6000
// #elif defined(TARGET_ARMV8A_ARM_CORTEX_A73)
// // common
// #define CACHE_LINE_SIZE 64
// // double
// #define D_PS 4  // panel size
// #define D_PLD 4  // 2 // GCD of panel length
// #define D_KC 320
// #define D_NC 256
// #define D_MC 6000
// // single
// #define S_PS 4
// #define S_PLD 4  // 2
// #define S_KC 256
// #define S_NC 128  // TODO these are just dummy
// #define S_MC 3000  // TODO these are just dummy
// #elif defined(TARGET_ARMV8A_ARM_CORTEX_A57)
// // common
// #define CACHE_LINE_SIZE 64
// // double
// #define D_PS 4  // panel size
// #define D_PLD 4  // 2 // GCD of panel length
// #define D_KC 128  // 224 //256 //192
// #define D_NC 72  // 40 //36 //48
// #define D_MC (4 * 192)  // 512 //488 //600
// // single
// #define S_PS 4
// #define S_PLD 4  // 2
// #define S_KC 256
// #define S_NC 128  // TODO these are just dummy
// #define S_MC 3000  // TODO these are just dummy
// #elif defined(TARGET_ARMV8A_ARM_CORTEX_A55)
// // common
// #define CACHE_LINE_SIZE 64
// // double
// #define D_PS 4  // panel size
// #define D_PLD 4  // 2 // GCD of panel length
// #define D_KC 224
// #define D_NC 160
// #define D_MC 6000
// // single
// #define S_PS 4
// #define S_PLD 4  // 2
// #define S_KC 256
// #define S_NC 128  // TODO these are just dummy
// #define S_MC 3000  // TODO these are just dummy
// #elif defined(TARGET_ARMV8A_ARM_CORTEX_A53)
// // common
// #define CACHE_LINE_SIZE 64
// // double
// #define D_PS 4  // panel size
// #define D_PLD 4  // 2 // GCD of panel length
// #define D_KC 160
// #define D_NC 128
// #define D_MC 6000
// // single
// #define S_PS 4
// #define S_PLD 4  // 2
// #define S_KC 256
// #define S_NC 128  // TODO these are just dummy
// #define S_MC 3000  // TODO these are just dummy
// #elif defined(TARGET_ARMV7A_ARM_CORTEX_A15)
// // common
// #define CACHE_LINE_SIZE 64
// // double
// #define D_PS 4  // panel size
// #define D_PLD 4  // 2 // GCD of panel length
// #define D_KC 256
// #define D_NC 128  // TODO these are just dummy
// #define D_MC 3000  // TODO these are just dummy
// // single
// #define S_PS 4
// #define S_PLD 4  // 2
// #define S_KC 256
// #define S_NC 128  // TODO these are just dummy
// #define S_MC 3000  // TODO these are just dummy
// #elif defined(TARGET_ARMV7A_ARM_CORTEX_A7)
// // common
// #define CACHE_LINE_SIZE 64
// // double
// #define D_PS 4  // panel size
// #define D_PLD 4  // 2 // GCD of panel length
// #define D_KC 256
// #define D_NC 128  // TODO these are just dummy
// #define D_MC 3000  // TODO these are just dummy
// // single
// #define S_PS 4
// #define S_PLD 4  // 2
// #define S_KC 256
// #define S_NC 128  // TODO these are just dummy
// #define S_MC 3000  // TODO these are just dummy
// #elif defined(TARGET_ARMV7A_ARM_CORTEX_A9)
// // common
// #define CACHE_LINE_SIZE 32
// // double
// #define D_PS 4  // panel size
// #define D_PLD 4  // 2 // GCD of panel length
// #define D_KC 256
// #define D_NC 128  // TODO these are just dummy
// #define D_MC 3000  // TODO these are just dummy
// // single
// #define S_PS 4
// #define S_PLD 4  // 2
// #define S_KC 256
// #define S_NC 128  // TODO these are just dummy
// #define S_MC 3000  // TODO these are just dummy


#else
#error "Unknown architecture"
#endif


#define D_CACHE_LINE_EL (CACHE_LINE_SIZE / D_EL_SIZE)
#define D_L1_CACHE_EL (L1_CACHE_SIZE / D_EL_SIZE)
#define D_L2_CACHE_EL (L2_CACHE_SIZE / D_EL_SIZE)
#define D_LLC_CACHE_EL (LLC_CACHE_SIZE / D_EL_SIZE)

#define S_CACHE_LINE_EL (CACHE_LINE_SIZE / S_EL_SIZE)
#define S_L1_CACHE_EL (L1_CACHE_SIZE / S_EL_SIZE)
#define S_L2_CACHE_EL (L2_CACHE_SIZE / S_EL_SIZE)
#define S_LLC_CACHE_EL (LLC_CACHE_SIZE / S_EL_SIZE)


#endif  // BLASFEO_BLOCK_SIZE_H_
