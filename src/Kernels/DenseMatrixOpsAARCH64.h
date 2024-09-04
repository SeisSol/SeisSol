
#include <Kernels/Precision.h>

#ifdef __aarch64__

// cf. https://stackoverflow.com/a/61248308

#if defined(DOUBLE_PRECISION)

#define DMO_INCREMENT 2
#define DMO_STREAM(IN, OUT) asm volatile ("ldnp x9,x10,[%[INaddr]]\n stnp x9,x10,[%[OUTaddr]]" :: [INaddr] "r" (IN), [OUTaddr] "r" (OUT) : "x9", "x10");

#elif defined(SINGLE_PRECISION)

#define DMO_INCREMENT 4
#define DMO_STREAM(IN, OUT) asm volatile ("ldnp x9,x10,[%[INaddr]]\n stnp x9,x10,[%[OUTaddr]]" :: [INaddr] "r" (IN), [OUTaddr] "r" (OUT) : "x9", "x10");

#else
#error no precision was defined
#endif

#endif
