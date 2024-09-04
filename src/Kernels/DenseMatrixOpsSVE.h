
#include <Kernels/Precision.h>

#ifdef __ARM_FEATURE_SVE
#include <arm_sve.h>

#if defined(DOUBLE_PRECISION)

#define DMO_INCREMENT svcntd()
#define DMO_STREAM(IN, OUT) svstnt1_f64(svptrue_b64(), OUT, svldnt1_f64(svptrue_b64(), IN));

#elif defined(SINGLE_PRECISION)

#define DMO_INCREMENT svcntw()
#define DMO_STREAM(IN, OUT) svstnt1_f32(svptrue_b32(), OUT, svldnt1_f32(svptrue_b32(), IN));

#else
#error no precision was defined
#endif

#endif
