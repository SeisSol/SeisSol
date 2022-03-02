#ifndef NDEBUG
extern long long libxsmm_num_total_flops;
extern long long pspamm_num_total_flops;
#endif
#if defined( __SSE3__) || defined(__MIC__)
#include <immintrin.h>
#endif
void libxsmm_m3_n1_k3_ldA3_ldB3_ldC3_alpha1_beta0_alignedA0_alignedC0_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $1, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $2, %%r10\n\t"
                       "vxorpd %%xmm15, %%xmm15, %%xmm15\n\t"
                       "vmovupd 0(%%rdi), %%xmm1\n\t"
                       "addq $24, %%rdi\n\t"
                       "vmovddup 0(%%rsi), %%xmm0\n\t"
                       "vfmadd231pd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovupd 0(%%rdi), %%xmm1\n\t"
                       "addq $24, %%rdi\n\t"
                       "vmovddup 8(%%rsi), %%xmm0\n\t"
                       "vfmadd231pd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovupd 0(%%rdi), %%xmm1\n\t"
                       "addq $24, %%rdi\n\t"
                       "vmovddup 16(%%rsi), %%xmm0\n\t"
                       "vfmadd231pd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovupd %%xmm15, 0(%%rdx)\n\t"
                       "addq $16, %%rdx\n\t"
                       "subq $56, %%rdi\n\t"
                       "cmpq $2, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $2, %%r10\n\t"
                       "34:\n\t"
                       "addq $1, %%r10\n\t"
                       "vxorpd %%xmm15, %%xmm15, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $24, %%rdi\n\t"
                       "vmovsd 0(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $24, %%rdi\n\t"
                       "vmovsd 8(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $24, %%rdi\n\t"
                       "vmovsd 16(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd %%xmm15, 0(%%rdx)\n\t"
                       "addq $8, %%rdx\n\t"
                       "subq $64, %%rdi\n\t"
                       "cmpq $3, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $0, %%rdx\n\t"
                       "addq $24, %%rsi\n\t"
                       "subq $24, %%rdi\n\t"
                       "cmpq $1, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 18;
#endif
}

void libxsmm_m3_n1_k3_ldA3_ldB3_ldC3_alpha1_beta1_alignedA0_alignedC0_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $1, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $2, %%r10\n\t"
                       "vmovupd 0(%%rdx), %%xmm15\n\t"
                       "vmovupd 0(%%rdi), %%xmm1\n\t"
                       "addq $24, %%rdi\n\t"
                       "vmovddup 0(%%rsi), %%xmm0\n\t"
                       "vfmadd231pd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovupd 0(%%rdi), %%xmm1\n\t"
                       "addq $24, %%rdi\n\t"
                       "vmovddup 8(%%rsi), %%xmm0\n\t"
                       "vfmadd231pd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovupd 0(%%rdi), %%xmm1\n\t"
                       "addq $24, %%rdi\n\t"
                       "vmovddup 16(%%rsi), %%xmm0\n\t"
                       "vfmadd231pd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovupd %%xmm15, 0(%%rdx)\n\t"
                       "addq $16, %%rdx\n\t"
                       "subq $56, %%rdi\n\t"
                       "cmpq $2, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $2, %%r10\n\t"
                       "34:\n\t"
                       "addq $1, %%r10\n\t"
                       "vmovsd 0(%%rdx), %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $24, %%rdi\n\t"
                       "vmovsd 0(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $24, %%rdi\n\t"
                       "vmovsd 8(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $24, %%rdi\n\t"
                       "vmovsd 16(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd %%xmm15, 0(%%rdx)\n\t"
                       "addq $8, %%rdx\n\t"
                       "subq $64, %%rdi\n\t"
                       "cmpq $3, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $0, %%rdx\n\t"
                       "addq $24, %%rsi\n\t"
                       "subq $24, %%rdi\n\t"
                       "cmpq $1, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 18;
#endif
}

void libxsmm_m56_n9_k343_ldA56_ldB344_ldC56_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vxorpd %%ymm4, %%ymm4, %%ymm4\n\t"
                       "vxorpd %%ymm5, %%ymm5, %%ymm5\n\t"
                       "vxorpd %%ymm6, %%ymm6, %%ymm6\n\t"
                       "vxorpd %%ymm7, %%ymm7, %%ymm7\n\t"
                       "vxorpd %%ymm8, %%ymm8, %%ymm8\n\t"
                       "vxorpd %%ymm9, %%ymm9, %%ymm9\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 2752(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 5504(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 2752(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 5504(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 2752(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 5504(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 2752(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 5504(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "cmpq $340, %%r12\n\t"
                       "jl 35b\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 2752(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 5504(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 2752(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 5504(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 2752(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 5504(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "subq $2744, %%rsi\n\t"
                       "vmovapd %%ymm4, 0(%%rdx)\n\t"
                       "vmovapd %%ymm5, 32(%%rdx)\n\t"
                       "vmovapd %%ymm6, 64(%%rdx)\n\t"
                       "vmovapd %%ymm7, 96(%%rdx)\n\t"
                       "vmovapd %%ymm8, 448(%%rdx)\n\t"
                       "vmovapd %%ymm9, 480(%%rdx)\n\t"
                       "vmovapd %%ymm10, 512(%%rdx)\n\t"
                       "vmovapd %%ymm11, 544(%%rdx)\n\t"
                       "vmovapd %%ymm12, 896(%%rdx)\n\t"
                       "vmovapd %%ymm13, 928(%%rdx)\n\t"
                       "vmovapd %%ymm14, 960(%%rdx)\n\t"
                       "vmovapd %%ymm15, 992(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $153536, %%rdi\n\t"
                       "cmpq $32, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $32, %%r10\n\t"
                       "34:\n\t"
                       "addq $12, %%r10\n\t"
                       "vxorpd %%ymm7, %%ymm7, %%ymm7\n\t"
                       "vxorpd %%ymm8, %%ymm8, %%ymm8\n\t"
                       "vxorpd %%ymm9, %%ymm9, %%ymm9\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 2752(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 5504(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 2752(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 5504(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 2752(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 5504(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 2752(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 5504(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "cmpq $340, %%r12\n\t"
                       "jl 35b\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 2752(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 5504(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 2752(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 5504(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 2752(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 5504(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "subq $2744, %%rsi\n\t"
                       "vmovapd %%ymm7, 0(%%rdx)\n\t"
                       "vmovapd %%ymm8, 32(%%rdx)\n\t"
                       "vmovapd %%ymm9, 64(%%rdx)\n\t"
                       "vmovapd %%ymm10, 448(%%rdx)\n\t"
                       "vmovapd %%ymm11, 480(%%rdx)\n\t"
                       "vmovapd %%ymm12, 512(%%rdx)\n\t"
                       "vmovapd %%ymm13, 896(%%rdx)\n\t"
                       "vmovapd %%ymm14, 928(%%rdx)\n\t"
                       "vmovapd %%ymm15, 960(%%rdx)\n\t"
                       "addq $96, %%rdx\n\t"
                       "subq $153568, %%rdi\n\t"
                       "cmpq $56, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $896, %%rdx\n\t"
                       "addq $8256, %%rsi\n\t"
                       "subq $448, %%rdi\n\t"
                       "cmpq $9, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 345744;
#endif
}

void libxsmm_m344_n9_k56_ldA344_ldB56_ldC344_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vxorpd %%ymm4, %%ymm4, %%ymm4\n\t"
                       "vxorpd %%ymm5, %%ymm5, %%ymm5\n\t"
                       "vxorpd %%ymm6, %%ymm6, %%ymm6\n\t"
                       "vxorpd %%ymm7, %%ymm7, %%ymm7\n\t"
                       "vxorpd %%ymm8, %%ymm8, %%ymm8\n\t"
                       "vxorpd %%ymm9, %%ymm9, %%ymm9\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $2752, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $2752, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $2752, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $2752, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "cmpq $56, %%r12\n\t"
                       "jl 35b\n\t"
                       "subq $448, %%rsi\n\t"
                       "vmovapd %%ymm4, 0(%%rdx)\n\t"
                       "vmovapd %%ymm5, 32(%%rdx)\n\t"
                       "vmovapd %%ymm6, 64(%%rdx)\n\t"
                       "vmovapd %%ymm7, 96(%%rdx)\n\t"
                       "vmovapd %%ymm8, 2752(%%rdx)\n\t"
                       "vmovapd %%ymm9, 2784(%%rdx)\n\t"
                       "vmovapd %%ymm10, 2816(%%rdx)\n\t"
                       "vmovapd %%ymm11, 2848(%%rdx)\n\t"
                       "vmovapd %%ymm12, 5504(%%rdx)\n\t"
                       "vmovapd %%ymm13, 5536(%%rdx)\n\t"
                       "vmovapd %%ymm14, 5568(%%rdx)\n\t"
                       "vmovapd %%ymm15, 5600(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $153984, %%rdi\n\t"
                       "cmpq $336, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $336, %%r10\n\t"
                       "34:\n\t"
                       "addq $8, %%r10\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $2752, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $2752, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $2752, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $2752, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "cmpq $56, %%r12\n\t"
                       "jl 35b\n\t"
                       "subq $448, %%rsi\n\t"
                       "vmovapd %%ymm10, 0(%%rdx)\n\t"
                       "vmovapd %%ymm11, 32(%%rdx)\n\t"
                       "vmovapd %%ymm12, 2752(%%rdx)\n\t"
                       "vmovapd %%ymm13, 2784(%%rdx)\n\t"
                       "vmovapd %%ymm14, 5504(%%rdx)\n\t"
                       "vmovapd %%ymm15, 5536(%%rdx)\n\t"
                       "addq $64, %%rdx\n\t"
                       "subq $154048, %%rdi\n\t"
                       "cmpq $344, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $5504, %%rdx\n\t"
                       "addq $1344, %%rsi\n\t"
                       "subq $2752, %%rdi\n\t"
                       "cmpq $9, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 346752;
#endif
}

void libxsmm_m36_n9_k9_ldA56_ldB9_ldC36_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vxorpd %%ymm4, %%ymm4, %%ymm4\n\t"
                       "vxorpd %%ymm5, %%ymm5, %%ymm5\n\t"
                       "vxorpd %%ymm6, %%ymm6, %%ymm6\n\t"
                       "vxorpd %%ymm7, %%ymm7, %%ymm7\n\t"
                       "vxorpd %%ymm8, %%ymm8, %%ymm8\n\t"
                       "vxorpd %%ymm9, %%ymm9, %%ymm9\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 168(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 176(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 184(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm4, 0(%%rdx)\n\t"
                       "vmovapd %%ymm5, 32(%%rdx)\n\t"
                       "vmovapd %%ymm6, 64(%%rdx)\n\t"
                       "vmovapd %%ymm7, 96(%%rdx)\n\t"
                       "vmovapd %%ymm8, 288(%%rdx)\n\t"
                       "vmovapd %%ymm9, 320(%%rdx)\n\t"
                       "vmovapd %%ymm10, 352(%%rdx)\n\t"
                       "vmovapd %%ymm11, 384(%%rdx)\n\t"
                       "vmovapd %%ymm12, 576(%%rdx)\n\t"
                       "vmovapd %%ymm13, 608(%%rdx)\n\t"
                       "vmovapd %%ymm14, 640(%%rdx)\n\t"
                       "vmovapd %%ymm15, 672(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $3904, %%rdi\n\t"
                       "cmpq $32, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $32, %%r10\n\t"
                       "34:\n\t"
                       "addq $4, %%r10\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 168(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 176(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 184(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm13, 0(%%rdx)\n\t"
                       "vmovapd %%ymm14, 288(%%rdx)\n\t"
                       "vmovapd %%ymm15, 576(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $4000, %%rdi\n\t"
                       "cmpq $36, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $576, %%rdx\n\t"
                       "addq $216, %%rsi\n\t"
                       "subq $288, %%rdi\n\t"
                       "cmpq $9, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 5832;
#endif
}

void libxsmm_m56_n9_k35_ldA56_ldB36_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vmovapd 0(%%rdx), %%ymm4\n\t"
                       "vmovapd 32(%%rdx), %%ymm5\n\t"
                       "vmovapd 64(%%rdx), %%ymm6\n\t"
                       "vmovapd 96(%%rdx), %%ymm7\n\t"
                       "vmovapd 448(%%rdx), %%ymm8\n\t"
                       "vmovapd 480(%%rdx), %%ymm9\n\t"
                       "vmovapd 512(%%rdx), %%ymm10\n\t"
                       "vmovapd 544(%%rdx), %%ymm11\n\t"
                       "vmovapd 896(%%rdx), %%ymm12\n\t"
                       "vmovapd 928(%%rdx), %%ymm13\n\t"
                       "vmovapd 960(%%rdx), %%ymm14\n\t"
                       "vmovapd 992(%%rdx), %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 576(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 576(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 576(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 576(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "cmpq $32, %%r12\n\t"
                       "jl 35b\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 576(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 576(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 576(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "subq $280, %%rsi\n\t"
                       "vmovapd %%ymm4, 0(%%rdx)\n\t"
                       "vmovapd %%ymm5, 32(%%rdx)\n\t"
                       "vmovapd %%ymm6, 64(%%rdx)\n\t"
                       "vmovapd %%ymm7, 96(%%rdx)\n\t"
                       "vmovapd %%ymm8, 448(%%rdx)\n\t"
                       "vmovapd %%ymm9, 480(%%rdx)\n\t"
                       "vmovapd %%ymm10, 512(%%rdx)\n\t"
                       "vmovapd %%ymm11, 544(%%rdx)\n\t"
                       "vmovapd %%ymm12, 896(%%rdx)\n\t"
                       "vmovapd %%ymm13, 928(%%rdx)\n\t"
                       "vmovapd %%ymm14, 960(%%rdx)\n\t"
                       "vmovapd %%ymm15, 992(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $15552, %%rdi\n\t"
                       "cmpq $32, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $32, %%r10\n\t"
                       "34:\n\t"
                       "addq $12, %%r10\n\t"
                       "vmovapd 0(%%rdx), %%ymm7\n\t"
                       "vmovapd 32(%%rdx), %%ymm8\n\t"
                       "vmovapd 64(%%rdx), %%ymm9\n\t"
                       "vmovapd 448(%%rdx), %%ymm10\n\t"
                       "vmovapd 480(%%rdx), %%ymm11\n\t"
                       "vmovapd 512(%%rdx), %%ymm12\n\t"
                       "vmovapd 896(%%rdx), %%ymm13\n\t"
                       "vmovapd 928(%%rdx), %%ymm14\n\t"
                       "vmovapd 960(%%rdx), %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 576(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 576(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 576(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 576(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "cmpq $32, %%r12\n\t"
                       "jl 35b\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 576(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 576(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 576(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "subq $280, %%rsi\n\t"
                       "vmovapd %%ymm7, 0(%%rdx)\n\t"
                       "vmovapd %%ymm8, 32(%%rdx)\n\t"
                       "vmovapd %%ymm9, 64(%%rdx)\n\t"
                       "vmovapd %%ymm10, 448(%%rdx)\n\t"
                       "vmovapd %%ymm11, 480(%%rdx)\n\t"
                       "vmovapd %%ymm12, 512(%%rdx)\n\t"
                       "vmovapd %%ymm13, 896(%%rdx)\n\t"
                       "vmovapd %%ymm14, 928(%%rdx)\n\t"
                       "vmovapd %%ymm15, 960(%%rdx)\n\t"
                       "addq $96, %%rdx\n\t"
                       "subq $15584, %%rdi\n\t"
                       "cmpq $56, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $896, %%rdx\n\t"
                       "addq $864, %%rsi\n\t"
                       "subq $448, %%rdi\n\t"
                       "cmpq $9, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 35280;
#endif
}

void libxsmm_m56_n6_k56_ldA56_ldB56_ldC56_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vxorpd %%ymm4, %%ymm4, %%ymm4\n\t"
                       "vxorpd %%ymm5, %%ymm5, %%ymm5\n\t"
                       "vxorpd %%ymm6, %%ymm6, %%ymm6\n\t"
                       "vxorpd %%ymm7, %%ymm7, %%ymm7\n\t"
                       "vxorpd %%ymm8, %%ymm8, %%ymm8\n\t"
                       "vxorpd %%ymm9, %%ymm9, %%ymm9\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "cmpq $56, %%r12\n\t"
                       "jl 35b\n\t"
                       "subq $448, %%rsi\n\t"
                       "vmovapd %%ymm4, 0(%%rdx)\n\t"
                       "vmovapd %%ymm5, 32(%%rdx)\n\t"
                       "vmovapd %%ymm6, 64(%%rdx)\n\t"
                       "vmovapd %%ymm7, 96(%%rdx)\n\t"
                       "vmovapd %%ymm8, 448(%%rdx)\n\t"
                       "vmovapd %%ymm9, 480(%%rdx)\n\t"
                       "vmovapd %%ymm10, 512(%%rdx)\n\t"
                       "vmovapd %%ymm11, 544(%%rdx)\n\t"
                       "vmovapd %%ymm12, 896(%%rdx)\n\t"
                       "vmovapd %%ymm13, 928(%%rdx)\n\t"
                       "vmovapd %%ymm14, 960(%%rdx)\n\t"
                       "vmovapd %%ymm15, 992(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $24960, %%rdi\n\t"
                       "cmpq $32, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $32, %%r10\n\t"
                       "34:\n\t"
                       "addq $12, %%r10\n\t"
                       "vxorpd %%ymm7, %%ymm7, %%ymm7\n\t"
                       "vxorpd %%ymm8, %%ymm8, %%ymm8\n\t"
                       "vxorpd %%ymm9, %%ymm9, %%ymm9\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "cmpq $56, %%r12\n\t"
                       "jl 35b\n\t"
                       "subq $448, %%rsi\n\t"
                       "vmovapd %%ymm7, 0(%%rdx)\n\t"
                       "vmovapd %%ymm8, 32(%%rdx)\n\t"
                       "vmovapd %%ymm9, 64(%%rdx)\n\t"
                       "vmovapd %%ymm10, 448(%%rdx)\n\t"
                       "vmovapd %%ymm11, 480(%%rdx)\n\t"
                       "vmovapd %%ymm12, 512(%%rdx)\n\t"
                       "vmovapd %%ymm13, 896(%%rdx)\n\t"
                       "vmovapd %%ymm14, 928(%%rdx)\n\t"
                       "vmovapd %%ymm15, 960(%%rdx)\n\t"
                       "addq $96, %%rdx\n\t"
                       "subq $24992, %%rdi\n\t"
                       "cmpq $56, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $896, %%rdx\n\t"
                       "addq $1344, %%rsi\n\t"
                       "subq $448, %%rdi\n\t"
                       "cmpq $6, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 37632;
#endif
}

void libxsmm_m56_n1_k56_ldA56_ldB56_ldC56_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $1, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm1\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm1\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm1\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm1\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm15\n\t"
                       "cmpq $56, %%r12\n\t"
                       "jl 35b\n\t"
                       "subq $448, %%rsi\n\t"
                       "vmovapd %%ymm12, 0(%%rdx)\n\t"
                       "vmovapd %%ymm13, 32(%%rdx)\n\t"
                       "vmovapd %%ymm14, 64(%%rdx)\n\t"
                       "vmovapd %%ymm15, 96(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $24960, %%rdi\n\t"
                       "cmpq $32, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $32, %%r10\n\t"
                       "34:\n\t"
                       "addq $12, %%r10\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm1\n\t"
                       "vmovapd 32(%%rdi), %%ymm2\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm13\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm1\n\t"
                       "vmovapd 32(%%rdi), %%ymm2\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm13\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm1\n\t"
                       "vmovapd 32(%%rdi), %%ymm2\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm13\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm1\n\t"
                       "vmovapd 32(%%rdi), %%ymm2\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm13\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm15\n\t"
                       "cmpq $56, %%r12\n\t"
                       "jl 35b\n\t"
                       "subq $448, %%rsi\n\t"
                       "vmovapd %%ymm13, 0(%%rdx)\n\t"
                       "vmovapd %%ymm14, 32(%%rdx)\n\t"
                       "vmovapd %%ymm15, 64(%%rdx)\n\t"
                       "addq $96, %%rdx\n\t"
                       "subq $24992, %%rdi\n\t"
                       "cmpq $56, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $0, %%rdx\n\t"
                       "addq $448, %%rsi\n\t"
                       "subq $448, %%rdi\n\t"
                       "cmpq $1, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 6272;
#endif
}

void libxsmm_m56_n1_k3_ldA56_ldB3_ldC56_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $1, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vmovapd 0(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm1\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vmovapd 0(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm1\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vmovapd 0(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm1\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm15\n\t"
                       "vmovapd %%ymm12, 0(%%rdx)\n\t"
                       "vmovapd %%ymm13, 32(%%rdx)\n\t"
                       "vmovapd %%ymm14, 64(%%rdx)\n\t"
                       "vmovapd %%ymm15, 96(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $1216, %%rdi\n\t"
                       "cmpq $32, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $32, %%r10\n\t"
                       "34:\n\t"
                       "addq $12, %%r10\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vmovapd 0(%%rdi), %%ymm1\n\t"
                       "vmovapd 32(%%rdi), %%ymm2\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm13\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vmovapd 0(%%rdi), %%ymm1\n\t"
                       "vmovapd 32(%%rdi), %%ymm2\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm13\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vmovapd 0(%%rdi), %%ymm1\n\t"
                       "vmovapd 32(%%rdi), %%ymm2\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm13\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm15\n\t"
                       "vmovapd %%ymm13, 0(%%rdx)\n\t"
                       "vmovapd %%ymm14, 32(%%rdx)\n\t"
                       "vmovapd %%ymm15, 64(%%rdx)\n\t"
                       "addq $96, %%rdx\n\t"
                       "subq $1248, %%rdi\n\t"
                       "cmpq $56, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $0, %%rdx\n\t"
                       "addq $24, %%rsi\n\t"
                       "subq $448, %%rdi\n\t"
                       "cmpq $1, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 336;
#endif
}

void libxsmm_m56_n1_k6_ldA56_ldB6_ldC56_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $1, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vmovapd 0(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm1\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vmovapd 0(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm1\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vmovapd 0(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm1\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vmovapd 0(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm1\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vmovapd 0(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm1\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vmovapd 0(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm1\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm15\n\t"
                       "vmovapd %%ymm12, 0(%%rdx)\n\t"
                       "vmovapd %%ymm13, 32(%%rdx)\n\t"
                       "vmovapd %%ymm14, 64(%%rdx)\n\t"
                       "vmovapd %%ymm15, 96(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $2560, %%rdi\n\t"
                       "cmpq $32, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $32, %%r10\n\t"
                       "34:\n\t"
                       "addq $12, %%r10\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vmovapd 0(%%rdi), %%ymm1\n\t"
                       "vmovapd 32(%%rdi), %%ymm2\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm13\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vmovapd 0(%%rdi), %%ymm1\n\t"
                       "vmovapd 32(%%rdi), %%ymm2\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm13\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vmovapd 0(%%rdi), %%ymm1\n\t"
                       "vmovapd 32(%%rdi), %%ymm2\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm13\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vmovapd 0(%%rdi), %%ymm1\n\t"
                       "vmovapd 32(%%rdi), %%ymm2\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm13\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vmovapd 0(%%rdi), %%ymm1\n\t"
                       "vmovapd 32(%%rdi), %%ymm2\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm13\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vmovapd 0(%%rdi), %%ymm1\n\t"
                       "vmovapd 32(%%rdi), %%ymm2\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm13\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm15\n\t"
                       "vmovapd %%ymm13, 0(%%rdx)\n\t"
                       "vmovapd %%ymm14, 32(%%rdx)\n\t"
                       "vmovapd %%ymm15, 64(%%rdx)\n\t"
                       "addq $96, %%rdx\n\t"
                       "subq $2592, %%rdi\n\t"
                       "cmpq $56, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $0, %%rdx\n\t"
                       "addq $48, %%rsi\n\t"
                       "subq $448, %%rdi\n\t"
                       "cmpq $1, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 672;
#endif
}

void libxsmm_m56_n6_k56_ldA56_ldB56_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vmovapd 0(%%rdx), %%ymm4\n\t"
                       "vmovapd 32(%%rdx), %%ymm5\n\t"
                       "vmovapd 64(%%rdx), %%ymm6\n\t"
                       "vmovapd 96(%%rdx), %%ymm7\n\t"
                       "vmovapd 448(%%rdx), %%ymm8\n\t"
                       "vmovapd 480(%%rdx), %%ymm9\n\t"
                       "vmovapd 512(%%rdx), %%ymm10\n\t"
                       "vmovapd 544(%%rdx), %%ymm11\n\t"
                       "vmovapd 896(%%rdx), %%ymm12\n\t"
                       "vmovapd 928(%%rdx), %%ymm13\n\t"
                       "vmovapd 960(%%rdx), %%ymm14\n\t"
                       "vmovapd 992(%%rdx), %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "cmpq $56, %%r12\n\t"
                       "jl 35b\n\t"
                       "subq $448, %%rsi\n\t"
                       "vmovapd %%ymm4, 0(%%rdx)\n\t"
                       "vmovapd %%ymm5, 32(%%rdx)\n\t"
                       "vmovapd %%ymm6, 64(%%rdx)\n\t"
                       "vmovapd %%ymm7, 96(%%rdx)\n\t"
                       "vmovapd %%ymm8, 448(%%rdx)\n\t"
                       "vmovapd %%ymm9, 480(%%rdx)\n\t"
                       "vmovapd %%ymm10, 512(%%rdx)\n\t"
                       "vmovapd %%ymm11, 544(%%rdx)\n\t"
                       "vmovapd %%ymm12, 896(%%rdx)\n\t"
                       "vmovapd %%ymm13, 928(%%rdx)\n\t"
                       "vmovapd %%ymm14, 960(%%rdx)\n\t"
                       "vmovapd %%ymm15, 992(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $24960, %%rdi\n\t"
                       "cmpq $32, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $32, %%r10\n\t"
                       "34:\n\t"
                       "addq $12, %%r10\n\t"
                       "vmovapd 0(%%rdx), %%ymm7\n\t"
                       "vmovapd 32(%%rdx), %%ymm8\n\t"
                       "vmovapd 64(%%rdx), %%ymm9\n\t"
                       "vmovapd 448(%%rdx), %%ymm10\n\t"
                       "vmovapd 480(%%rdx), %%ymm11\n\t"
                       "vmovapd 512(%%rdx), %%ymm12\n\t"
                       "vmovapd 896(%%rdx), %%ymm13\n\t"
                       "vmovapd 928(%%rdx), %%ymm14\n\t"
                       "vmovapd 960(%%rdx), %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "cmpq $56, %%r12\n\t"
                       "jl 35b\n\t"
                       "subq $448, %%rsi\n\t"
                       "vmovapd %%ymm7, 0(%%rdx)\n\t"
                       "vmovapd %%ymm8, 32(%%rdx)\n\t"
                       "vmovapd %%ymm9, 64(%%rdx)\n\t"
                       "vmovapd %%ymm10, 448(%%rdx)\n\t"
                       "vmovapd %%ymm11, 480(%%rdx)\n\t"
                       "vmovapd %%ymm12, 512(%%rdx)\n\t"
                       "vmovapd %%ymm13, 896(%%rdx)\n\t"
                       "vmovapd %%ymm14, 928(%%rdx)\n\t"
                       "vmovapd %%ymm15, 960(%%rdx)\n\t"
                       "addq $96, %%rdx\n\t"
                       "subq $24992, %%rdi\n\t"
                       "cmpq $56, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $896, %%rdx\n\t"
                       "addq $1344, %%rsi\n\t"
                       "subq $448, %%rdi\n\t"
                       "cmpq $6, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 37632;
#endif
}

void libxsmm_m9_n1_k9_ldA9_ldB9_ldC9_alpha1_beta1_alignedA0_alignedC0_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $1, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $8, %%r10\n\t"
                       "vmovupd 0(%%rdx), %%ymm14\n\t"
                       "vmovupd 32(%%rdx), %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "vmovupd 32(%%rdi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "addq $72, %%rdi\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "vmovupd 32(%%rdi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "addq $72, %%rdi\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "vmovupd 32(%%rdi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "addq $72, %%rdi\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "vmovupd 32(%%rdi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "addq $72, %%rdi\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "vmovupd 32(%%rdi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "addq $72, %%rdi\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "vmovupd 32(%%rdi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "addq $72, %%rdi\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "vmovupd 32(%%rdi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "addq $72, %%rdi\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "vmovupd 32(%%rdi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "addq $72, %%rdi\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "vmovupd 32(%%rdi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "addq $72, %%rdi\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm15\n\t"
                       "vmovupd %%ymm14, 0(%%rdx)\n\t"
                       "vmovupd %%ymm15, 32(%%rdx)\n\t"
                       "addq $64, %%rdx\n\t"
                       "subq $584, %%rdi\n\t"
                       "cmpq $8, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $8, %%r10\n\t"
                       "34:\n\t"
                       "addq $1, %%r10\n\t"
                       "vmovsd 0(%%rdx), %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $72, %%rdi\n\t"
                       "vmovsd 0(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $72, %%rdi\n\t"
                       "vmovsd 8(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $72, %%rdi\n\t"
                       "vmovsd 16(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $72, %%rdi\n\t"
                       "vmovsd 24(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $72, %%rdi\n\t"
                       "vmovsd 32(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $72, %%rdi\n\t"
                       "vmovsd 40(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $72, %%rdi\n\t"
                       "vmovsd 48(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $72, %%rdi\n\t"
                       "vmovsd 56(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $72, %%rdi\n\t"
                       "vmovsd 64(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd %%xmm15, 0(%%rdx)\n\t"
                       "addq $8, %%rdx\n\t"
                       "subq $640, %%rdi\n\t"
                       "cmpq $9, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $0, %%rdx\n\t"
                       "addq $72, %%rsi\n\t"
                       "subq $72, %%rdi\n\t"
                       "cmpq $1, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 162;
#endif
}

void libxsmm_m9_n1_k9_ldA9_ldB9_ldC9_alpha1_beta0_alignedA0_alignedC0_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $1, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $8, %%r10\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "vmovupd 32(%%rdi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "addq $72, %%rdi\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "vmovupd 32(%%rdi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "addq $72, %%rdi\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "vmovupd 32(%%rdi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "addq $72, %%rdi\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "vmovupd 32(%%rdi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "addq $72, %%rdi\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "vmovupd 32(%%rdi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "addq $72, %%rdi\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "vmovupd 32(%%rdi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "addq $72, %%rdi\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "vmovupd 32(%%rdi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "addq $72, %%rdi\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "vmovupd 32(%%rdi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "addq $72, %%rdi\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "vmovupd 32(%%rdi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "addq $72, %%rdi\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm15\n\t"
                       "vmovupd %%ymm14, 0(%%rdx)\n\t"
                       "vmovupd %%ymm15, 32(%%rdx)\n\t"
                       "addq $64, %%rdx\n\t"
                       "subq $584, %%rdi\n\t"
                       "cmpq $8, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $8, %%r10\n\t"
                       "34:\n\t"
                       "addq $1, %%r10\n\t"
                       "vxorpd %%xmm15, %%xmm15, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $72, %%rdi\n\t"
                       "vmovsd 0(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $72, %%rdi\n\t"
                       "vmovsd 8(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $72, %%rdi\n\t"
                       "vmovsd 16(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $72, %%rdi\n\t"
                       "vmovsd 24(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $72, %%rdi\n\t"
                       "vmovsd 32(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $72, %%rdi\n\t"
                       "vmovsd 40(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $72, %%rdi\n\t"
                       "vmovsd 48(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $72, %%rdi\n\t"
                       "vmovsd 56(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $72, %%rdi\n\t"
                       "vmovsd 64(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd %%xmm15, 0(%%rdx)\n\t"
                       "addq $8, %%rdx\n\t"
                       "subq $640, %%rdi\n\t"
                       "cmpq $9, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $0, %%rdx\n\t"
                       "addq $72, %%rsi\n\t"
                       "subq $72, %%rdi\n\t"
                       "cmpq $1, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 162;
#endif
}

void libxsmm_m6_n1_k9_ldA6_ldB9_ldC6_alpha1_beta0_alignedA0_alignedC0_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $1, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $4, %%r10\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "addq $48, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "addq $48, %%rdi\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "addq $48, %%rdi\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "addq $48, %%rdi\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "addq $48, %%rdi\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "addq $48, %%rdi\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "addq $48, %%rdi\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "addq $48, %%rdi\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "addq $48, %%rdi\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm15\n\t"
                       "vmovupd %%ymm15, 0(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $400, %%rdi\n\t"
                       "cmpq $4, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $4, %%r10\n\t"
                       "34:\n\t"
                       "addq $2, %%r10\n\t"
                       "vxorpd %%xmm15, %%xmm15, %%xmm15\n\t"
                       "vmovupd 0(%%rdi), %%xmm1\n\t"
                       "addq $48, %%rdi\n\t"
                       "vmovddup 0(%%rsi), %%xmm0\n\t"
                       "vfmadd231pd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovupd 0(%%rdi), %%xmm1\n\t"
                       "addq $48, %%rdi\n\t"
                       "vmovddup 8(%%rsi), %%xmm0\n\t"
                       "vfmadd231pd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovupd 0(%%rdi), %%xmm1\n\t"
                       "addq $48, %%rdi\n\t"
                       "vmovddup 16(%%rsi), %%xmm0\n\t"
                       "vfmadd231pd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovupd 0(%%rdi), %%xmm1\n\t"
                       "addq $48, %%rdi\n\t"
                       "vmovddup 24(%%rsi), %%xmm0\n\t"
                       "vfmadd231pd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovupd 0(%%rdi), %%xmm1\n\t"
                       "addq $48, %%rdi\n\t"
                       "vmovddup 32(%%rsi), %%xmm0\n\t"
                       "vfmadd231pd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovupd 0(%%rdi), %%xmm1\n\t"
                       "addq $48, %%rdi\n\t"
                       "vmovddup 40(%%rsi), %%xmm0\n\t"
                       "vfmadd231pd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovupd 0(%%rdi), %%xmm1\n\t"
                       "addq $48, %%rdi\n\t"
                       "vmovddup 48(%%rsi), %%xmm0\n\t"
                       "vfmadd231pd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovupd 0(%%rdi), %%xmm1\n\t"
                       "addq $48, %%rdi\n\t"
                       "vmovddup 56(%%rsi), %%xmm0\n\t"
                       "vfmadd231pd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovupd 0(%%rdi), %%xmm1\n\t"
                       "addq $48, %%rdi\n\t"
                       "vmovddup 64(%%rsi), %%xmm0\n\t"
                       "vfmadd231pd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovupd %%xmm15, 0(%%rdx)\n\t"
                       "addq $16, %%rdx\n\t"
                       "subq $416, %%rdi\n\t"
                       "cmpq $6, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $0, %%rdx\n\t"
                       "addq $72, %%rsi\n\t"
                       "subq $48, %%rdi\n\t"
                       "cmpq $1, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 108;
#endif
}

void libxsmm_m9_n1_k6_ldA9_ldB6_ldC9_alpha1_beta0_alignedA0_alignedC0_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $1, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $8, %%r10\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "vmovupd 32(%%rdi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "addq $72, %%rdi\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "vmovupd 32(%%rdi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "addq $72, %%rdi\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "vmovupd 32(%%rdi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "addq $72, %%rdi\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "vmovupd 32(%%rdi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "addq $72, %%rdi\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "vmovupd 32(%%rdi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "addq $72, %%rdi\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vmovupd 0(%%rdi), %%ymm1\n\t"
                       "vmovupd 32(%%rdi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm1, %%ymm0, %%ymm14\n\t"
                       "addq $72, %%rdi\n\t"
                       "vfmadd231pd %%ymm2, %%ymm0, %%ymm15\n\t"
                       "vmovupd %%ymm14, 0(%%rdx)\n\t"
                       "vmovupd %%ymm15, 32(%%rdx)\n\t"
                       "addq $64, %%rdx\n\t"
                       "subq $368, %%rdi\n\t"
                       "cmpq $8, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $8, %%r10\n\t"
                       "34:\n\t"
                       "addq $1, %%r10\n\t"
                       "vxorpd %%xmm15, %%xmm15, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $72, %%rdi\n\t"
                       "vmovsd 0(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $72, %%rdi\n\t"
                       "vmovsd 8(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $72, %%rdi\n\t"
                       "vmovsd 16(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $72, %%rdi\n\t"
                       "vmovsd 24(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $72, %%rdi\n\t"
                       "vmovsd 32(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm1\n\t"
                       "addq $72, %%rdi\n\t"
                       "vmovsd 40(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm1, %%xmm0, %%xmm15\n\t"
                       "vmovsd %%xmm15, 0(%%rdx)\n\t"
                       "addq $8, %%rdx\n\t"
                       "subq $424, %%rdi\n\t"
                       "cmpq $9, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $0, %%rdx\n\t"
                       "addq $48, %%rsi\n\t"
                       "subq $72, %%rdi\n\t"
                       "cmpq $1, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 108;
#endif
}

void libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vxorpd %%ymm4, %%ymm4, %%ymm4\n\t"
                       "vxorpd %%ymm5, %%ymm5, %%ymm5\n\t"
                       "vxorpd %%ymm6, %%ymm6, %%ymm6\n\t"
                       "vxorpd %%ymm7, %%ymm7, %%ymm7\n\t"
                       "vxorpd %%ymm8, %%ymm8, %%ymm8\n\t"
                       "vxorpd %%ymm9, %%ymm9, %%ymm9\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "cmpq $56, %%r12\n\t"
                       "jl 35b\n\t"
                       "subq $448, %%rsi\n\t"
                       "vmovapd %%ymm4, 0(%%rdx)\n\t"
                       "vmovapd %%ymm5, 32(%%rdx)\n\t"
                       "vmovapd %%ymm6, 64(%%rdx)\n\t"
                       "vmovapd %%ymm7, 96(%%rdx)\n\t"
                       "vmovapd %%ymm8, 192(%%rdx)\n\t"
                       "vmovapd %%ymm9, 224(%%rdx)\n\t"
                       "vmovapd %%ymm10, 256(%%rdx)\n\t"
                       "vmovapd %%ymm11, 288(%%rdx)\n\t"
                       "vmovapd %%ymm12, 384(%%rdx)\n\t"
                       "vmovapd %%ymm13, 416(%%rdx)\n\t"
                       "vmovapd %%ymm14, 448(%%rdx)\n\t"
                       "vmovapd %%ymm15, 480(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $10624, %%rdi\n\t"
                       "cmpq $16, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $16, %%r10\n\t"
                       "34:\n\t"
                       "addq $8, %%r10\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "cmpq $56, %%r12\n\t"
                       "jl 35b\n\t"
                       "subq $448, %%rsi\n\t"
                       "vmovapd %%ymm10, 0(%%rdx)\n\t"
                       "vmovapd %%ymm11, 32(%%rdx)\n\t"
                       "vmovapd %%ymm12, 192(%%rdx)\n\t"
                       "vmovapd %%ymm13, 224(%%rdx)\n\t"
                       "vmovapd %%ymm14, 384(%%rdx)\n\t"
                       "vmovapd %%ymm15, 416(%%rdx)\n\t"
                       "addq $64, %%rdx\n\t"
                       "subq $10688, %%rdi\n\t"
                       "cmpq $24, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $384, %%rdx\n\t"
                       "addq $1344, %%rsi\n\t"
                       "subq $192, %%rdi\n\t"
                       "cmpq $9, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 24192;
#endif
}

void libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vxorpd %%ymm4, %%ymm4, %%ymm4\n\t"
                       "vxorpd %%ymm5, %%ymm5, %%ymm5\n\t"
                       "vxorpd %%ymm6, %%ymm6, %%ymm6\n\t"
                       "vxorpd %%ymm7, %%ymm7, %%ymm7\n\t"
                       "vxorpd %%ymm8, %%ymm8, %%ymm8\n\t"
                       "vxorpd %%ymm9, %%ymm9, %%ymm9\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 168(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 176(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 184(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm4, 0(%%rdx)\n\t"
                       "vmovapd %%ymm5, 32(%%rdx)\n\t"
                       "vmovapd %%ymm6, 64(%%rdx)\n\t"
                       "vmovapd %%ymm7, 96(%%rdx)\n\t"
                       "vmovapd %%ymm8, 192(%%rdx)\n\t"
                       "vmovapd %%ymm9, 224(%%rdx)\n\t"
                       "vmovapd %%ymm10, 256(%%rdx)\n\t"
                       "vmovapd %%ymm11, 288(%%rdx)\n\t"
                       "vmovapd %%ymm12, 384(%%rdx)\n\t"
                       "vmovapd %%ymm13, 416(%%rdx)\n\t"
                       "vmovapd %%ymm14, 448(%%rdx)\n\t"
                       "vmovapd %%ymm15, 480(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $1600, %%rdi\n\t"
                       "cmpq $16, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $16, %%r10\n\t"
                       "34:\n\t"
                       "addq $8, %%r10\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 168(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 176(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 184(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm10, 0(%%rdx)\n\t"
                       "vmovapd %%ymm11, 32(%%rdx)\n\t"
                       "vmovapd %%ymm12, 192(%%rdx)\n\t"
                       "vmovapd %%ymm13, 224(%%rdx)\n\t"
                       "vmovapd %%ymm14, 384(%%rdx)\n\t"
                       "vmovapd %%ymm15, 416(%%rdx)\n\t"
                       "addq $64, %%rdx\n\t"
                       "subq $1664, %%rdi\n\t"
                       "cmpq $24, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $384, %%rdx\n\t"
                       "addq $216, %%rsi\n\t"
                       "subq $192, %%rdi\n\t"
                       "cmpq $9, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 3888;
#endif
}

void libxsmmsparse_38193e2bc92c2e551f139883f9580e10_m56_n9_k21_ldA0_ldB24_ldC56_alpha1_beta1_alignedA0_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_n = 0;
  #pragma nounroll_and_jam
  for ( l_n = 0; l_n < 9; l_n++) {
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b0 = _mm256_broadcast_sd(&B[(l_n*24)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b0 = _mm_loaddup_pd(&B[(l_n*24)+0]);
#endif
    __m128d c0_0 = _mm_load_sd(&C[(l_n*56)+0]);
    __m128d a0_0 = _mm_load_sd(&A[0]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, b0));
#endif
    _mm_store_sd(&C[(l_n*56)+0], c0_0);
    __m128d c0_1 = _mm_load_sd(&C[(l_n*56)+3]);
    __m128d a0_1 = _mm_load_sd(&A[1]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, b0));
#endif
    _mm_store_sd(&C[(l_n*56)+3], c0_1);
    __m128d c0_2 = _mm_load_sd(&C[(l_n*56)+9]);
    __m128d a0_2 = _mm_load_sd(&A[2]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, b0));
#endif
    _mm_store_sd(&C[(l_n*56)+9], c0_2);
    __m128d c0_3 = _mm_load_sd(&C[(l_n*56)+19]);
    __m128d a0_3 = _mm_load_sd(&A[3]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, b0));
#endif
    _mm_store_sd(&C[(l_n*56)+19], c0_3);
    __m128d c0_4 = _mm_load_sd(&C[(l_n*56)+34]);
    __m128d a0_4 = _mm_load_sd(&A[4]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, b0));
#endif
    _mm_store_sd(&C[(l_n*56)+34], c0_4);
    __m128d c0_5 = _mm_load_sd(&C[(l_n*56)+55]);
    __m128d a0_5 = _mm_load_sd(&A[5]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_5 = _mm_add_sd(c0_5, _mm_mul_sd(a0_5, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_5 = _mm_add_sd(c0_5, _mm_mul_sd(a0_5, b0));
#endif
    _mm_store_sd(&C[(l_n*56)+55], c0_5);
#else
    C[(l_n*56)+0] += A[0] * B[(l_n*24)+0];
    C[(l_n*56)+3] += A[1] * B[(l_n*24)+0];
    C[(l_n*56)+9] += A[2] * B[(l_n*24)+0];
    C[(l_n*56)+19] += A[3] * B[(l_n*24)+0];
    C[(l_n*56)+34] += A[4] * B[(l_n*24)+0];
    C[(l_n*56)+55] += A[5] * B[(l_n*24)+0];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b1 = _mm256_broadcast_sd(&B[(l_n*24)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b1 = _mm_loaddup_pd(&B[(l_n*24)+1]);
#endif
    __m128d c1_0 = _mm_loadu_pd(&C[(l_n*56)+1]);
    __m128d a1_0 = _mm_loadu_pd(&A[6]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_0 = _mm_add_pd(c1_0, _mm_mul_pd(a1_0, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_0 = _mm_add_pd(c1_0, _mm_mul_pd(a1_0, b1));
#endif
    _mm_storeu_pd(&C[(l_n*56)+1], c1_0);
    __m128d c1_2 = _mm_loadu_pd(&C[(l_n*56)+7]);
    __m128d a1_2 = _mm_loadu_pd(&A[8]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_2 = _mm_add_pd(c1_2, _mm_mul_pd(a1_2, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_2 = _mm_add_pd(c1_2, _mm_mul_pd(a1_2, b1));
#endif
    _mm_storeu_pd(&C[(l_n*56)+7], c1_2);
    __m128d c1_4 = _mm_loadu_pd(&C[(l_n*56)+17]);
    __m128d a1_4 = _mm_loadu_pd(&A[10]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_4 = _mm_add_pd(c1_4, _mm_mul_pd(a1_4, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_4 = _mm_add_pd(c1_4, _mm_mul_pd(a1_4, b1));
#endif
    _mm_storeu_pd(&C[(l_n*56)+17], c1_4);
    __m128d c1_6 = _mm_loadu_pd(&C[(l_n*56)+32]);
    __m128d a1_6 = _mm_loadu_pd(&A[12]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_6 = _mm_add_pd(c1_6, _mm_mul_pd(a1_6, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_6 = _mm_add_pd(c1_6, _mm_mul_pd(a1_6, b1));
#endif
    _mm_storeu_pd(&C[(l_n*56)+32], c1_6);
    __m128d c1_8 = _mm_loadu_pd(&C[(l_n*56)+53]);
    __m128d a1_8 = _mm_loadu_pd(&A[14]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_8 = _mm_add_pd(c1_8, _mm_mul_pd(a1_8, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_8 = _mm_add_pd(c1_8, _mm_mul_pd(a1_8, b1));
#endif
    _mm_storeu_pd(&C[(l_n*56)+53], c1_8);
#else
    C[(l_n*56)+1] += A[6] * B[(l_n*24)+1];
    C[(l_n*56)+2] += A[7] * B[(l_n*24)+1];
    C[(l_n*56)+7] += A[8] * B[(l_n*24)+1];
    C[(l_n*56)+8] += A[9] * B[(l_n*24)+1];
    C[(l_n*56)+17] += A[10] * B[(l_n*24)+1];
    C[(l_n*56)+18] += A[11] * B[(l_n*24)+1];
    C[(l_n*56)+32] += A[12] * B[(l_n*24)+1];
    C[(l_n*56)+33] += A[13] * B[(l_n*24)+1];
    C[(l_n*56)+53] += A[14] * B[(l_n*24)+1];
    C[(l_n*56)+54] += A[15] * B[(l_n*24)+1];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b2 = _mm256_broadcast_sd(&B[(l_n*24)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b2 = _mm_loaddup_pd(&B[(l_n*24)+2]);
#endif
    __m128d c2_0 = _mm_loadu_pd(&C[(l_n*56)+1]);
    __m128d a2_0 = _mm_loadu_pd(&A[16]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_0 = _mm_add_pd(c2_0, _mm_mul_pd(a2_0, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_0 = _mm_add_pd(c2_0, _mm_mul_pd(a2_0, b2));
#endif
    _mm_storeu_pd(&C[(l_n*56)+1], c2_0);
    __m128d c2_2 = _mm_loadu_pd(&C[(l_n*56)+7]);
    __m128d a2_2 = _mm_loadu_pd(&A[18]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_2 = _mm_add_pd(c2_2, _mm_mul_pd(a2_2, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_2 = _mm_add_pd(c2_2, _mm_mul_pd(a2_2, b2));
#endif
    _mm_storeu_pd(&C[(l_n*56)+7], c2_2);
    __m128d c2_4 = _mm_loadu_pd(&C[(l_n*56)+17]);
    __m128d a2_4 = _mm_loadu_pd(&A[20]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_4 = _mm_add_pd(c2_4, _mm_mul_pd(a2_4, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_4 = _mm_add_pd(c2_4, _mm_mul_pd(a2_4, b2));
#endif
    _mm_storeu_pd(&C[(l_n*56)+17], c2_4);
    __m128d c2_6 = _mm_loadu_pd(&C[(l_n*56)+32]);
    __m128d a2_6 = _mm_loadu_pd(&A[22]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_6 = _mm_add_pd(c2_6, _mm_mul_pd(a2_6, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_6 = _mm_add_pd(c2_6, _mm_mul_pd(a2_6, b2));
#endif
    _mm_storeu_pd(&C[(l_n*56)+32], c2_6);
    __m128d c2_8 = _mm_loadu_pd(&C[(l_n*56)+53]);
    __m128d a2_8 = _mm_loadu_pd(&A[24]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_8 = _mm_add_pd(c2_8, _mm_mul_pd(a2_8, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_8 = _mm_add_pd(c2_8, _mm_mul_pd(a2_8, b2));
#endif
    _mm_storeu_pd(&C[(l_n*56)+53], c2_8);
#else
    C[(l_n*56)+1] += A[16] * B[(l_n*24)+2];
    C[(l_n*56)+2] += A[17] * B[(l_n*24)+2];
    C[(l_n*56)+7] += A[18] * B[(l_n*24)+2];
    C[(l_n*56)+8] += A[19] * B[(l_n*24)+2];
    C[(l_n*56)+17] += A[20] * B[(l_n*24)+2];
    C[(l_n*56)+18] += A[21] * B[(l_n*24)+2];
    C[(l_n*56)+32] += A[22] * B[(l_n*24)+2];
    C[(l_n*56)+33] += A[23] * B[(l_n*24)+2];
    C[(l_n*56)+53] += A[24] * B[(l_n*24)+2];
    C[(l_n*56)+54] += A[25] * B[(l_n*24)+2];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b3 = _mm256_broadcast_sd(&B[(l_n*24)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b3 = _mm_loaddup_pd(&B[(l_n*24)+3]);
#endif
    __m128d c3_0 = _mm_loadu_pd(&C[(l_n*56)+4]);
    __m128d a3_0 = _mm_loadu_pd(&A[26]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_0 = _mm_add_pd(c3_0, _mm_mul_pd(a3_0, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_0 = _mm_add_pd(c3_0, _mm_mul_pd(a3_0, b3));
#endif
    _mm_storeu_pd(&C[(l_n*56)+4], c3_0);
    __m128d c3_2 = _mm_load_sd(&C[(l_n*56)+6]);
    __m128d a3_2 = _mm_load_sd(&A[28]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, b3));
#endif
    _mm_store_sd(&C[(l_n*56)+6], c3_2);
    __m128d c3_3 = _mm_loadu_pd(&C[(l_n*56)+14]);
    __m128d a3_3 = _mm_loadu_pd(&A[29]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_3 = _mm_add_pd(c3_3, _mm_mul_pd(a3_3, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_3 = _mm_add_pd(c3_3, _mm_mul_pd(a3_3, b3));
#endif
    _mm_storeu_pd(&C[(l_n*56)+14], c3_3);
    __m128d c3_5 = _mm_load_sd(&C[(l_n*56)+16]);
    __m128d a3_5 = _mm_load_sd(&A[31]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_5 = _mm_add_sd(c3_5, _mm_mul_sd(a3_5, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_5 = _mm_add_sd(c3_5, _mm_mul_sd(a3_5, b3));
#endif
    _mm_store_sd(&C[(l_n*56)+16], c3_5);
    __m128d c3_6 = _mm_loadu_pd(&C[(l_n*56)+29]);
    __m128d a3_6 = _mm_loadu_pd(&A[32]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_6 = _mm_add_pd(c3_6, _mm_mul_pd(a3_6, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_6 = _mm_add_pd(c3_6, _mm_mul_pd(a3_6, b3));
#endif
    _mm_storeu_pd(&C[(l_n*56)+29], c3_6);
    __m128d c3_8 = _mm_load_sd(&C[(l_n*56)+31]);
    __m128d a3_8 = _mm_load_sd(&A[34]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_8 = _mm_add_sd(c3_8, _mm_mul_sd(a3_8, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_8 = _mm_add_sd(c3_8, _mm_mul_sd(a3_8, b3));
#endif
    _mm_store_sd(&C[(l_n*56)+31], c3_8);
    __m128d c3_9 = _mm_loadu_pd(&C[(l_n*56)+50]);
    __m128d a3_9 = _mm_loadu_pd(&A[35]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_9 = _mm_add_pd(c3_9, _mm_mul_pd(a3_9, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_9 = _mm_add_pd(c3_9, _mm_mul_pd(a3_9, b3));
#endif
    _mm_storeu_pd(&C[(l_n*56)+50], c3_9);
    __m128d c3_11 = _mm_load_sd(&C[(l_n*56)+52]);
    __m128d a3_11 = _mm_load_sd(&A[37]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_11 = _mm_add_sd(c3_11, _mm_mul_sd(a3_11, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_11 = _mm_add_sd(c3_11, _mm_mul_sd(a3_11, b3));
#endif
    _mm_store_sd(&C[(l_n*56)+52], c3_11);
#else
    C[(l_n*56)+4] += A[26] * B[(l_n*24)+3];
    C[(l_n*56)+5] += A[27] * B[(l_n*24)+3];
    C[(l_n*56)+6] += A[28] * B[(l_n*24)+3];
    C[(l_n*56)+14] += A[29] * B[(l_n*24)+3];
    C[(l_n*56)+15] += A[30] * B[(l_n*24)+3];
    C[(l_n*56)+16] += A[31] * B[(l_n*24)+3];
    C[(l_n*56)+29] += A[32] * B[(l_n*24)+3];
    C[(l_n*56)+30] += A[33] * B[(l_n*24)+3];
    C[(l_n*56)+31] += A[34] * B[(l_n*24)+3];
    C[(l_n*56)+50] += A[35] * B[(l_n*24)+3];
    C[(l_n*56)+51] += A[36] * B[(l_n*24)+3];
    C[(l_n*56)+52] += A[37] * B[(l_n*24)+3];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b4 = _mm256_broadcast_sd(&B[(l_n*24)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b4 = _mm_loaddup_pd(&B[(l_n*24)+4]);
#endif
    __m128d c4_0 = _mm_loadu_pd(&C[(l_n*56)+4]);
    __m128d a4_0 = _mm_loadu_pd(&A[38]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_0 = _mm_add_pd(c4_0, _mm_mul_pd(a4_0, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_0 = _mm_add_pd(c4_0, _mm_mul_pd(a4_0, b4));
#endif
    _mm_storeu_pd(&C[(l_n*56)+4], c4_0);
    __m128d c4_2 = _mm_load_sd(&C[(l_n*56)+6]);
    __m128d a4_2 = _mm_load_sd(&A[40]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, b4));
#endif
    _mm_store_sd(&C[(l_n*56)+6], c4_2);
    __m128d c4_3 = _mm_loadu_pd(&C[(l_n*56)+14]);
    __m128d a4_3 = _mm_loadu_pd(&A[41]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_3 = _mm_add_pd(c4_3, _mm_mul_pd(a4_3, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_3 = _mm_add_pd(c4_3, _mm_mul_pd(a4_3, b4));
#endif
    _mm_storeu_pd(&C[(l_n*56)+14], c4_3);
    __m128d c4_5 = _mm_load_sd(&C[(l_n*56)+16]);
    __m128d a4_5 = _mm_load_sd(&A[43]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_5 = _mm_add_sd(c4_5, _mm_mul_sd(a4_5, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_5 = _mm_add_sd(c4_5, _mm_mul_sd(a4_5, b4));
#endif
    _mm_store_sd(&C[(l_n*56)+16], c4_5);
    __m128d c4_6 = _mm_loadu_pd(&C[(l_n*56)+29]);
    __m128d a4_6 = _mm_loadu_pd(&A[44]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_6 = _mm_add_pd(c4_6, _mm_mul_pd(a4_6, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_6 = _mm_add_pd(c4_6, _mm_mul_pd(a4_6, b4));
#endif
    _mm_storeu_pd(&C[(l_n*56)+29], c4_6);
    __m128d c4_8 = _mm_load_sd(&C[(l_n*56)+31]);
    __m128d a4_8 = _mm_load_sd(&A[46]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_8 = _mm_add_sd(c4_8, _mm_mul_sd(a4_8, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_8 = _mm_add_sd(c4_8, _mm_mul_sd(a4_8, b4));
#endif
    _mm_store_sd(&C[(l_n*56)+31], c4_8);
    __m128d c4_9 = _mm_loadu_pd(&C[(l_n*56)+50]);
    __m128d a4_9 = _mm_loadu_pd(&A[47]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_9 = _mm_add_pd(c4_9, _mm_mul_pd(a4_9, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_9 = _mm_add_pd(c4_9, _mm_mul_pd(a4_9, b4));
#endif
    _mm_storeu_pd(&C[(l_n*56)+50], c4_9);
    __m128d c4_11 = _mm_load_sd(&C[(l_n*56)+52]);
    __m128d a4_11 = _mm_load_sd(&A[49]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_11 = _mm_add_sd(c4_11, _mm_mul_sd(a4_11, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_11 = _mm_add_sd(c4_11, _mm_mul_sd(a4_11, b4));
#endif
    _mm_store_sd(&C[(l_n*56)+52], c4_11);
#else
    C[(l_n*56)+4] += A[38] * B[(l_n*24)+4];
    C[(l_n*56)+5] += A[39] * B[(l_n*24)+4];
    C[(l_n*56)+6] += A[40] * B[(l_n*24)+4];
    C[(l_n*56)+14] += A[41] * B[(l_n*24)+4];
    C[(l_n*56)+15] += A[42] * B[(l_n*24)+4];
    C[(l_n*56)+16] += A[43] * B[(l_n*24)+4];
    C[(l_n*56)+29] += A[44] * B[(l_n*24)+4];
    C[(l_n*56)+30] += A[45] * B[(l_n*24)+4];
    C[(l_n*56)+31] += A[46] * B[(l_n*24)+4];
    C[(l_n*56)+50] += A[47] * B[(l_n*24)+4];
    C[(l_n*56)+51] += A[48] * B[(l_n*24)+4];
    C[(l_n*56)+52] += A[49] * B[(l_n*24)+4];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b5 = _mm256_broadcast_sd(&B[(l_n*24)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b5 = _mm_loaddup_pd(&B[(l_n*24)+5]);
#endif
    __m128d c5_0 = _mm_loadu_pd(&C[(l_n*56)+4]);
    __m128d a5_0 = _mm_loadu_pd(&A[50]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_0 = _mm_add_pd(c5_0, _mm_mul_pd(a5_0, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_0 = _mm_add_pd(c5_0, _mm_mul_pd(a5_0, b5));
#endif
    _mm_storeu_pd(&C[(l_n*56)+4], c5_0);
    __m128d c5_2 = _mm_load_sd(&C[(l_n*56)+6]);
    __m128d a5_2 = _mm_load_sd(&A[52]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, b5));
#endif
    _mm_store_sd(&C[(l_n*56)+6], c5_2);
    __m128d c5_3 = _mm_loadu_pd(&C[(l_n*56)+14]);
    __m128d a5_3 = _mm_loadu_pd(&A[53]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_3 = _mm_add_pd(c5_3, _mm_mul_pd(a5_3, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_3 = _mm_add_pd(c5_3, _mm_mul_pd(a5_3, b5));
#endif
    _mm_storeu_pd(&C[(l_n*56)+14], c5_3);
    __m128d c5_5 = _mm_load_sd(&C[(l_n*56)+16]);
    __m128d a5_5 = _mm_load_sd(&A[55]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_5 = _mm_add_sd(c5_5, _mm_mul_sd(a5_5, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_5 = _mm_add_sd(c5_5, _mm_mul_sd(a5_5, b5));
#endif
    _mm_store_sd(&C[(l_n*56)+16], c5_5);
    __m128d c5_6 = _mm_loadu_pd(&C[(l_n*56)+29]);
    __m128d a5_6 = _mm_loadu_pd(&A[56]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_6 = _mm_add_pd(c5_6, _mm_mul_pd(a5_6, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_6 = _mm_add_pd(c5_6, _mm_mul_pd(a5_6, b5));
#endif
    _mm_storeu_pd(&C[(l_n*56)+29], c5_6);
    __m128d c5_8 = _mm_load_sd(&C[(l_n*56)+31]);
    __m128d a5_8 = _mm_load_sd(&A[58]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_8 = _mm_add_sd(c5_8, _mm_mul_sd(a5_8, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_8 = _mm_add_sd(c5_8, _mm_mul_sd(a5_8, b5));
#endif
    _mm_store_sd(&C[(l_n*56)+31], c5_8);
    __m128d c5_9 = _mm_loadu_pd(&C[(l_n*56)+50]);
    __m128d a5_9 = _mm_loadu_pd(&A[59]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_9 = _mm_add_pd(c5_9, _mm_mul_pd(a5_9, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_9 = _mm_add_pd(c5_9, _mm_mul_pd(a5_9, b5));
#endif
    _mm_storeu_pd(&C[(l_n*56)+50], c5_9);
    __m128d c5_11 = _mm_load_sd(&C[(l_n*56)+52]);
    __m128d a5_11 = _mm_load_sd(&A[61]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_11 = _mm_add_sd(c5_11, _mm_mul_sd(a5_11, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_11 = _mm_add_sd(c5_11, _mm_mul_sd(a5_11, b5));
#endif
    _mm_store_sd(&C[(l_n*56)+52], c5_11);
#else
    C[(l_n*56)+4] += A[50] * B[(l_n*24)+5];
    C[(l_n*56)+5] += A[51] * B[(l_n*24)+5];
    C[(l_n*56)+6] += A[52] * B[(l_n*24)+5];
    C[(l_n*56)+14] += A[53] * B[(l_n*24)+5];
    C[(l_n*56)+15] += A[54] * B[(l_n*24)+5];
    C[(l_n*56)+16] += A[55] * B[(l_n*24)+5];
    C[(l_n*56)+29] += A[56] * B[(l_n*24)+5];
    C[(l_n*56)+30] += A[57] * B[(l_n*24)+5];
    C[(l_n*56)+31] += A[58] * B[(l_n*24)+5];
    C[(l_n*56)+50] += A[59] * B[(l_n*24)+5];
    C[(l_n*56)+51] += A[60] * B[(l_n*24)+5];
    C[(l_n*56)+52] += A[61] * B[(l_n*24)+5];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b6 = _mm256_broadcast_sd(&B[(l_n*24)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b6 = _mm_loaddup_pd(&B[(l_n*24)+6]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c6_0 = _mm256_loadu_pd(&C[(l_n*56)+10]);
    __m256d a6_0 = _mm256_loadu_pd(&A[62]);
    c6_0 = _mm256_add_pd(c6_0, _mm256_mul_pd(a6_0, b6));
    _mm256_storeu_pd(&C[(l_n*56)+10], c6_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c6_0 = _mm_loadu_pd(&C[(l_n*56)+10]);
    __m128d a6_0 = _mm_loadu_pd(&A[62]);
    c6_0 = _mm_add_pd(c6_0, _mm_mul_pd(a6_0, b6));
    _mm_storeu_pd(&C[(l_n*56)+10], c6_0);
    __m128d c6_2 = _mm_loadu_pd(&C[(l_n*56)+12]);
    __m128d a6_2 = _mm_loadu_pd(&A[64]);
    c6_2 = _mm_add_pd(c6_2, _mm_mul_pd(a6_2, b6));
    _mm_storeu_pd(&C[(l_n*56)+12], c6_2);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c6_4 = _mm256_loadu_pd(&C[(l_n*56)+25]);
    __m256d a6_4 = _mm256_loadu_pd(&A[66]);
    c6_4 = _mm256_add_pd(c6_4, _mm256_mul_pd(a6_4, b6));
    _mm256_storeu_pd(&C[(l_n*56)+25], c6_4);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c6_4 = _mm_loadu_pd(&C[(l_n*56)+25]);
    __m128d a6_4 = _mm_loadu_pd(&A[66]);
    c6_4 = _mm_add_pd(c6_4, _mm_mul_pd(a6_4, b6));
    _mm_storeu_pd(&C[(l_n*56)+25], c6_4);
    __m128d c6_6 = _mm_loadu_pd(&C[(l_n*56)+27]);
    __m128d a6_6 = _mm_loadu_pd(&A[68]);
    c6_6 = _mm_add_pd(c6_6, _mm_mul_pd(a6_6, b6));
    _mm_storeu_pd(&C[(l_n*56)+27], c6_6);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c6_8 = _mm256_loadu_pd(&C[(l_n*56)+46]);
    __m256d a6_8 = _mm256_loadu_pd(&A[70]);
    c6_8 = _mm256_add_pd(c6_8, _mm256_mul_pd(a6_8, b6));
    _mm256_storeu_pd(&C[(l_n*56)+46], c6_8);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c6_8 = _mm_loadu_pd(&C[(l_n*56)+46]);
    __m128d a6_8 = _mm_loadu_pd(&A[70]);
    c6_8 = _mm_add_pd(c6_8, _mm_mul_pd(a6_8, b6));
    _mm_storeu_pd(&C[(l_n*56)+46], c6_8);
    __m128d c6_10 = _mm_loadu_pd(&C[(l_n*56)+48]);
    __m128d a6_10 = _mm_loadu_pd(&A[72]);
    c6_10 = _mm_add_pd(c6_10, _mm_mul_pd(a6_10, b6));
    _mm_storeu_pd(&C[(l_n*56)+48], c6_10);
#endif
#else
    C[(l_n*56)+10] += A[62] * B[(l_n*24)+6];
    C[(l_n*56)+11] += A[63] * B[(l_n*24)+6];
    C[(l_n*56)+12] += A[64] * B[(l_n*24)+6];
    C[(l_n*56)+13] += A[65] * B[(l_n*24)+6];
    C[(l_n*56)+25] += A[66] * B[(l_n*24)+6];
    C[(l_n*56)+26] += A[67] * B[(l_n*24)+6];
    C[(l_n*56)+27] += A[68] * B[(l_n*24)+6];
    C[(l_n*56)+28] += A[69] * B[(l_n*24)+6];
    C[(l_n*56)+46] += A[70] * B[(l_n*24)+6];
    C[(l_n*56)+47] += A[71] * B[(l_n*24)+6];
    C[(l_n*56)+48] += A[72] * B[(l_n*24)+6];
    C[(l_n*56)+49] += A[73] * B[(l_n*24)+6];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b7 = _mm256_broadcast_sd(&B[(l_n*24)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b7 = _mm_loaddup_pd(&B[(l_n*24)+7]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c7_0 = _mm256_loadu_pd(&C[(l_n*56)+10]);
    __m256d a7_0 = _mm256_loadu_pd(&A[74]);
    c7_0 = _mm256_add_pd(c7_0, _mm256_mul_pd(a7_0, b7));
    _mm256_storeu_pd(&C[(l_n*56)+10], c7_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c7_0 = _mm_loadu_pd(&C[(l_n*56)+10]);
    __m128d a7_0 = _mm_loadu_pd(&A[74]);
    c7_0 = _mm_add_pd(c7_0, _mm_mul_pd(a7_0, b7));
    _mm_storeu_pd(&C[(l_n*56)+10], c7_0);
    __m128d c7_2 = _mm_loadu_pd(&C[(l_n*56)+12]);
    __m128d a7_2 = _mm_loadu_pd(&A[76]);
    c7_2 = _mm_add_pd(c7_2, _mm_mul_pd(a7_2, b7));
    _mm_storeu_pd(&C[(l_n*56)+12], c7_2);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c7_4 = _mm256_loadu_pd(&C[(l_n*56)+25]);
    __m256d a7_4 = _mm256_loadu_pd(&A[78]);
    c7_4 = _mm256_add_pd(c7_4, _mm256_mul_pd(a7_4, b7));
    _mm256_storeu_pd(&C[(l_n*56)+25], c7_4);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c7_4 = _mm_loadu_pd(&C[(l_n*56)+25]);
    __m128d a7_4 = _mm_loadu_pd(&A[78]);
    c7_4 = _mm_add_pd(c7_4, _mm_mul_pd(a7_4, b7));
    _mm_storeu_pd(&C[(l_n*56)+25], c7_4);
    __m128d c7_6 = _mm_loadu_pd(&C[(l_n*56)+27]);
    __m128d a7_6 = _mm_loadu_pd(&A[80]);
    c7_6 = _mm_add_pd(c7_6, _mm_mul_pd(a7_6, b7));
    _mm_storeu_pd(&C[(l_n*56)+27], c7_6);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c7_8 = _mm256_loadu_pd(&C[(l_n*56)+46]);
    __m256d a7_8 = _mm256_loadu_pd(&A[82]);
    c7_8 = _mm256_add_pd(c7_8, _mm256_mul_pd(a7_8, b7));
    _mm256_storeu_pd(&C[(l_n*56)+46], c7_8);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c7_8 = _mm_loadu_pd(&C[(l_n*56)+46]);
    __m128d a7_8 = _mm_loadu_pd(&A[82]);
    c7_8 = _mm_add_pd(c7_8, _mm_mul_pd(a7_8, b7));
    _mm_storeu_pd(&C[(l_n*56)+46], c7_8);
    __m128d c7_10 = _mm_loadu_pd(&C[(l_n*56)+48]);
    __m128d a7_10 = _mm_loadu_pd(&A[84]);
    c7_10 = _mm_add_pd(c7_10, _mm_mul_pd(a7_10, b7));
    _mm_storeu_pd(&C[(l_n*56)+48], c7_10);
#endif
#else
    C[(l_n*56)+10] += A[74] * B[(l_n*24)+7];
    C[(l_n*56)+11] += A[75] * B[(l_n*24)+7];
    C[(l_n*56)+12] += A[76] * B[(l_n*24)+7];
    C[(l_n*56)+13] += A[77] * B[(l_n*24)+7];
    C[(l_n*56)+25] += A[78] * B[(l_n*24)+7];
    C[(l_n*56)+26] += A[79] * B[(l_n*24)+7];
    C[(l_n*56)+27] += A[80] * B[(l_n*24)+7];
    C[(l_n*56)+28] += A[81] * B[(l_n*24)+7];
    C[(l_n*56)+46] += A[82] * B[(l_n*24)+7];
    C[(l_n*56)+47] += A[83] * B[(l_n*24)+7];
    C[(l_n*56)+48] += A[84] * B[(l_n*24)+7];
    C[(l_n*56)+49] += A[85] * B[(l_n*24)+7];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b8 = _mm256_broadcast_sd(&B[(l_n*24)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b8 = _mm_loaddup_pd(&B[(l_n*24)+8]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c8_0 = _mm256_loadu_pd(&C[(l_n*56)+10]);
    __m256d a8_0 = _mm256_loadu_pd(&A[86]);
    c8_0 = _mm256_add_pd(c8_0, _mm256_mul_pd(a8_0, b8));
    _mm256_storeu_pd(&C[(l_n*56)+10], c8_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c8_0 = _mm_loadu_pd(&C[(l_n*56)+10]);
    __m128d a8_0 = _mm_loadu_pd(&A[86]);
    c8_0 = _mm_add_pd(c8_0, _mm_mul_pd(a8_0, b8));
    _mm_storeu_pd(&C[(l_n*56)+10], c8_0);
    __m128d c8_2 = _mm_loadu_pd(&C[(l_n*56)+12]);
    __m128d a8_2 = _mm_loadu_pd(&A[88]);
    c8_2 = _mm_add_pd(c8_2, _mm_mul_pd(a8_2, b8));
    _mm_storeu_pd(&C[(l_n*56)+12], c8_2);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c8_4 = _mm256_loadu_pd(&C[(l_n*56)+25]);
    __m256d a8_4 = _mm256_loadu_pd(&A[90]);
    c8_4 = _mm256_add_pd(c8_4, _mm256_mul_pd(a8_4, b8));
    _mm256_storeu_pd(&C[(l_n*56)+25], c8_4);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c8_4 = _mm_loadu_pd(&C[(l_n*56)+25]);
    __m128d a8_4 = _mm_loadu_pd(&A[90]);
    c8_4 = _mm_add_pd(c8_4, _mm_mul_pd(a8_4, b8));
    _mm_storeu_pd(&C[(l_n*56)+25], c8_4);
    __m128d c8_6 = _mm_loadu_pd(&C[(l_n*56)+27]);
    __m128d a8_6 = _mm_loadu_pd(&A[92]);
    c8_6 = _mm_add_pd(c8_6, _mm_mul_pd(a8_6, b8));
    _mm_storeu_pd(&C[(l_n*56)+27], c8_6);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c8_8 = _mm256_loadu_pd(&C[(l_n*56)+46]);
    __m256d a8_8 = _mm256_loadu_pd(&A[94]);
    c8_8 = _mm256_add_pd(c8_8, _mm256_mul_pd(a8_8, b8));
    _mm256_storeu_pd(&C[(l_n*56)+46], c8_8);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c8_8 = _mm_loadu_pd(&C[(l_n*56)+46]);
    __m128d a8_8 = _mm_loadu_pd(&A[94]);
    c8_8 = _mm_add_pd(c8_8, _mm_mul_pd(a8_8, b8));
    _mm_storeu_pd(&C[(l_n*56)+46], c8_8);
    __m128d c8_10 = _mm_loadu_pd(&C[(l_n*56)+48]);
    __m128d a8_10 = _mm_loadu_pd(&A[96]);
    c8_10 = _mm_add_pd(c8_10, _mm_mul_pd(a8_10, b8));
    _mm_storeu_pd(&C[(l_n*56)+48], c8_10);
#endif
#else
    C[(l_n*56)+10] += A[86] * B[(l_n*24)+8];
    C[(l_n*56)+11] += A[87] * B[(l_n*24)+8];
    C[(l_n*56)+12] += A[88] * B[(l_n*24)+8];
    C[(l_n*56)+13] += A[89] * B[(l_n*24)+8];
    C[(l_n*56)+25] += A[90] * B[(l_n*24)+8];
    C[(l_n*56)+26] += A[91] * B[(l_n*24)+8];
    C[(l_n*56)+27] += A[92] * B[(l_n*24)+8];
    C[(l_n*56)+28] += A[93] * B[(l_n*24)+8];
    C[(l_n*56)+46] += A[94] * B[(l_n*24)+8];
    C[(l_n*56)+47] += A[95] * B[(l_n*24)+8];
    C[(l_n*56)+48] += A[96] * B[(l_n*24)+8];
    C[(l_n*56)+49] += A[97] * B[(l_n*24)+8];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b9 = _mm256_broadcast_sd(&B[(l_n*24)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b9 = _mm_loaddup_pd(&B[(l_n*24)+9]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c9_0 = _mm256_loadu_pd(&C[(l_n*56)+10]);
    __m256d a9_0 = _mm256_loadu_pd(&A[98]);
    c9_0 = _mm256_add_pd(c9_0, _mm256_mul_pd(a9_0, b9));
    _mm256_storeu_pd(&C[(l_n*56)+10], c9_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c9_0 = _mm_loadu_pd(&C[(l_n*56)+10]);
    __m128d a9_0 = _mm_loadu_pd(&A[98]);
    c9_0 = _mm_add_pd(c9_0, _mm_mul_pd(a9_0, b9));
    _mm_storeu_pd(&C[(l_n*56)+10], c9_0);
    __m128d c9_2 = _mm_loadu_pd(&C[(l_n*56)+12]);
    __m128d a9_2 = _mm_loadu_pd(&A[100]);
    c9_2 = _mm_add_pd(c9_2, _mm_mul_pd(a9_2, b9));
    _mm_storeu_pd(&C[(l_n*56)+12], c9_2);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c9_4 = _mm256_loadu_pd(&C[(l_n*56)+25]);
    __m256d a9_4 = _mm256_loadu_pd(&A[102]);
    c9_4 = _mm256_add_pd(c9_4, _mm256_mul_pd(a9_4, b9));
    _mm256_storeu_pd(&C[(l_n*56)+25], c9_4);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c9_4 = _mm_loadu_pd(&C[(l_n*56)+25]);
    __m128d a9_4 = _mm_loadu_pd(&A[102]);
    c9_4 = _mm_add_pd(c9_4, _mm_mul_pd(a9_4, b9));
    _mm_storeu_pd(&C[(l_n*56)+25], c9_4);
    __m128d c9_6 = _mm_loadu_pd(&C[(l_n*56)+27]);
    __m128d a9_6 = _mm_loadu_pd(&A[104]);
    c9_6 = _mm_add_pd(c9_6, _mm_mul_pd(a9_6, b9));
    _mm_storeu_pd(&C[(l_n*56)+27], c9_6);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c9_8 = _mm256_loadu_pd(&C[(l_n*56)+46]);
    __m256d a9_8 = _mm256_loadu_pd(&A[106]);
    c9_8 = _mm256_add_pd(c9_8, _mm256_mul_pd(a9_8, b9));
    _mm256_storeu_pd(&C[(l_n*56)+46], c9_8);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c9_8 = _mm_loadu_pd(&C[(l_n*56)+46]);
    __m128d a9_8 = _mm_loadu_pd(&A[106]);
    c9_8 = _mm_add_pd(c9_8, _mm_mul_pd(a9_8, b9));
    _mm_storeu_pd(&C[(l_n*56)+46], c9_8);
    __m128d c9_10 = _mm_loadu_pd(&C[(l_n*56)+48]);
    __m128d a9_10 = _mm_loadu_pd(&A[108]);
    c9_10 = _mm_add_pd(c9_10, _mm_mul_pd(a9_10, b9));
    _mm_storeu_pd(&C[(l_n*56)+48], c9_10);
#endif
#else
    C[(l_n*56)+10] += A[98] * B[(l_n*24)+9];
    C[(l_n*56)+11] += A[99] * B[(l_n*24)+9];
    C[(l_n*56)+12] += A[100] * B[(l_n*24)+9];
    C[(l_n*56)+13] += A[101] * B[(l_n*24)+9];
    C[(l_n*56)+25] += A[102] * B[(l_n*24)+9];
    C[(l_n*56)+26] += A[103] * B[(l_n*24)+9];
    C[(l_n*56)+27] += A[104] * B[(l_n*24)+9];
    C[(l_n*56)+28] += A[105] * B[(l_n*24)+9];
    C[(l_n*56)+46] += A[106] * B[(l_n*24)+9];
    C[(l_n*56)+47] += A[107] * B[(l_n*24)+9];
    C[(l_n*56)+48] += A[108] * B[(l_n*24)+9];
    C[(l_n*56)+49] += A[109] * B[(l_n*24)+9];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b10 = _mm256_broadcast_sd(&B[(l_n*24)+10]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b10 = _mm_loaddup_pd(&B[(l_n*24)+10]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c10_0 = _mm256_loadu_pd(&C[(l_n*56)+20]);
    __m256d a10_0 = _mm256_loadu_pd(&A[110]);
    c10_0 = _mm256_add_pd(c10_0, _mm256_mul_pd(a10_0, b10));
    _mm256_storeu_pd(&C[(l_n*56)+20], c10_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c10_0 = _mm_loadu_pd(&C[(l_n*56)+20]);
    __m128d a10_0 = _mm_loadu_pd(&A[110]);
    c10_0 = _mm_add_pd(c10_0, _mm_mul_pd(a10_0, b10));
    _mm_storeu_pd(&C[(l_n*56)+20], c10_0);
    __m128d c10_2 = _mm_loadu_pd(&C[(l_n*56)+22]);
    __m128d a10_2 = _mm_loadu_pd(&A[112]);
    c10_2 = _mm_add_pd(c10_2, _mm_mul_pd(a10_2, b10));
    _mm_storeu_pd(&C[(l_n*56)+22], c10_2);
#endif
    __m128d c10_4 = _mm_load_sd(&C[(l_n*56)+24]);
    __m128d a10_4 = _mm_load_sd(&A[114]);
#if defined(__SSE3__) && defined(__AVX__)
    c10_4 = _mm_add_sd(c10_4, _mm_mul_sd(a10_4, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c10_4 = _mm_add_sd(c10_4, _mm_mul_sd(a10_4, b10));
#endif
    _mm_store_sd(&C[(l_n*56)+24], c10_4);
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c10_5 = _mm256_loadu_pd(&C[(l_n*56)+41]);
    __m256d a10_5 = _mm256_loadu_pd(&A[115]);
    c10_5 = _mm256_add_pd(c10_5, _mm256_mul_pd(a10_5, b10));
    _mm256_storeu_pd(&C[(l_n*56)+41], c10_5);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c10_5 = _mm_loadu_pd(&C[(l_n*56)+41]);
    __m128d a10_5 = _mm_loadu_pd(&A[115]);
    c10_5 = _mm_add_pd(c10_5, _mm_mul_pd(a10_5, b10));
    _mm_storeu_pd(&C[(l_n*56)+41], c10_5);
    __m128d c10_7 = _mm_loadu_pd(&C[(l_n*56)+43]);
    __m128d a10_7 = _mm_loadu_pd(&A[117]);
    c10_7 = _mm_add_pd(c10_7, _mm_mul_pd(a10_7, b10));
    _mm_storeu_pd(&C[(l_n*56)+43], c10_7);
#endif
    __m128d c10_9 = _mm_load_sd(&C[(l_n*56)+45]);
    __m128d a10_9 = _mm_load_sd(&A[119]);
#if defined(__SSE3__) && defined(__AVX__)
    c10_9 = _mm_add_sd(c10_9, _mm_mul_sd(a10_9, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c10_9 = _mm_add_sd(c10_9, _mm_mul_sd(a10_9, b10));
#endif
    _mm_store_sd(&C[(l_n*56)+45], c10_9);
#else
    C[(l_n*56)+20] += A[110] * B[(l_n*24)+10];
    C[(l_n*56)+21] += A[111] * B[(l_n*24)+10];
    C[(l_n*56)+22] += A[112] * B[(l_n*24)+10];
    C[(l_n*56)+23] += A[113] * B[(l_n*24)+10];
    C[(l_n*56)+24] += A[114] * B[(l_n*24)+10];
    C[(l_n*56)+41] += A[115] * B[(l_n*24)+10];
    C[(l_n*56)+42] += A[116] * B[(l_n*24)+10];
    C[(l_n*56)+43] += A[117] * B[(l_n*24)+10];
    C[(l_n*56)+44] += A[118] * B[(l_n*24)+10];
    C[(l_n*56)+45] += A[119] * B[(l_n*24)+10];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b11 = _mm256_broadcast_sd(&B[(l_n*24)+11]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b11 = _mm_loaddup_pd(&B[(l_n*24)+11]);
#endif
    __m128d c11_0 = _mm_loadu_pd(&C[(l_n*56)+20]);
    __m128d a11_0 = _mm_loadu_pd(&A[120]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_0 = _mm_add_pd(c11_0, _mm_mul_pd(a11_0, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_0 = _mm_add_pd(c11_0, _mm_mul_pd(a11_0, b11));
#endif
    _mm_storeu_pd(&C[(l_n*56)+20], c11_0);
    __m128d c11_2 = _mm_load_sd(&C[(l_n*56)+22]);
    __m128d a11_2 = _mm_load_sd(&A[122]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_2 = _mm_add_sd(c11_2, _mm_mul_sd(a11_2, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_2 = _mm_add_sd(c11_2, _mm_mul_sd(a11_2, b11));
#endif
    _mm_store_sd(&C[(l_n*56)+22], c11_2);
    __m128d c11_3 = _mm_load_sd(&C[(l_n*56)+24]);
    __m128d a11_3 = _mm_load_sd(&A[123]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_3 = _mm_add_sd(c11_3, _mm_mul_sd(a11_3, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_3 = _mm_add_sd(c11_3, _mm_mul_sd(a11_3, b11));
#endif
    _mm_store_sd(&C[(l_n*56)+24], c11_3);
    __m128d c11_4 = _mm_loadu_pd(&C[(l_n*56)+41]);
    __m128d a11_4 = _mm_loadu_pd(&A[124]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_4 = _mm_add_pd(c11_4, _mm_mul_pd(a11_4, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_4 = _mm_add_pd(c11_4, _mm_mul_pd(a11_4, b11));
#endif
    _mm_storeu_pd(&C[(l_n*56)+41], c11_4);
    __m128d c11_6 = _mm_load_sd(&C[(l_n*56)+43]);
    __m128d a11_6 = _mm_load_sd(&A[126]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_6 = _mm_add_sd(c11_6, _mm_mul_sd(a11_6, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_6 = _mm_add_sd(c11_6, _mm_mul_sd(a11_6, b11));
#endif
    _mm_store_sd(&C[(l_n*56)+43], c11_6);
    __m128d c11_7 = _mm_load_sd(&C[(l_n*56)+45]);
    __m128d a11_7 = _mm_load_sd(&A[127]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_7 = _mm_add_sd(c11_7, _mm_mul_sd(a11_7, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_7 = _mm_add_sd(c11_7, _mm_mul_sd(a11_7, b11));
#endif
    _mm_store_sd(&C[(l_n*56)+45], c11_7);
#else
    C[(l_n*56)+20] += A[120] * B[(l_n*24)+11];
    C[(l_n*56)+21] += A[121] * B[(l_n*24)+11];
    C[(l_n*56)+22] += A[122] * B[(l_n*24)+11];
    C[(l_n*56)+24] += A[123] * B[(l_n*24)+11];
    C[(l_n*56)+41] += A[124] * B[(l_n*24)+11];
    C[(l_n*56)+42] += A[125] * B[(l_n*24)+11];
    C[(l_n*56)+43] += A[126] * B[(l_n*24)+11];
    C[(l_n*56)+45] += A[127] * B[(l_n*24)+11];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b12 = _mm256_broadcast_sd(&B[(l_n*24)+12]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b12 = _mm_loaddup_pd(&B[(l_n*24)+12]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c12_0 = _mm256_loadu_pd(&C[(l_n*56)+20]);
    __m256d a12_0 = _mm256_loadu_pd(&A[128]);
    c12_0 = _mm256_add_pd(c12_0, _mm256_mul_pd(a12_0, b12));
    _mm256_storeu_pd(&C[(l_n*56)+20], c12_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c12_0 = _mm_loadu_pd(&C[(l_n*56)+20]);
    __m128d a12_0 = _mm_loadu_pd(&A[128]);
    c12_0 = _mm_add_pd(c12_0, _mm_mul_pd(a12_0, b12));
    _mm_storeu_pd(&C[(l_n*56)+20], c12_0);
    __m128d c12_2 = _mm_loadu_pd(&C[(l_n*56)+22]);
    __m128d a12_2 = _mm_loadu_pd(&A[130]);
    c12_2 = _mm_add_pd(c12_2, _mm_mul_pd(a12_2, b12));
    _mm_storeu_pd(&C[(l_n*56)+22], c12_2);
#endif
    __m128d c12_4 = _mm_load_sd(&C[(l_n*56)+24]);
    __m128d a12_4 = _mm_load_sd(&A[132]);
#if defined(__SSE3__) && defined(__AVX__)
    c12_4 = _mm_add_sd(c12_4, _mm_mul_sd(a12_4, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c12_4 = _mm_add_sd(c12_4, _mm_mul_sd(a12_4, b12));
#endif
    _mm_store_sd(&C[(l_n*56)+24], c12_4);
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c12_5 = _mm256_loadu_pd(&C[(l_n*56)+41]);
    __m256d a12_5 = _mm256_loadu_pd(&A[133]);
    c12_5 = _mm256_add_pd(c12_5, _mm256_mul_pd(a12_5, b12));
    _mm256_storeu_pd(&C[(l_n*56)+41], c12_5);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c12_5 = _mm_loadu_pd(&C[(l_n*56)+41]);
    __m128d a12_5 = _mm_loadu_pd(&A[133]);
    c12_5 = _mm_add_pd(c12_5, _mm_mul_pd(a12_5, b12));
    _mm_storeu_pd(&C[(l_n*56)+41], c12_5);
    __m128d c12_7 = _mm_loadu_pd(&C[(l_n*56)+43]);
    __m128d a12_7 = _mm_loadu_pd(&A[135]);
    c12_7 = _mm_add_pd(c12_7, _mm_mul_pd(a12_7, b12));
    _mm_storeu_pd(&C[(l_n*56)+43], c12_7);
#endif
    __m128d c12_9 = _mm_load_sd(&C[(l_n*56)+45]);
    __m128d a12_9 = _mm_load_sd(&A[137]);
#if defined(__SSE3__) && defined(__AVX__)
    c12_9 = _mm_add_sd(c12_9, _mm_mul_sd(a12_9, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c12_9 = _mm_add_sd(c12_9, _mm_mul_sd(a12_9, b12));
#endif
    _mm_store_sd(&C[(l_n*56)+45], c12_9);
#else
    C[(l_n*56)+20] += A[128] * B[(l_n*24)+12];
    C[(l_n*56)+21] += A[129] * B[(l_n*24)+12];
    C[(l_n*56)+22] += A[130] * B[(l_n*24)+12];
    C[(l_n*56)+23] += A[131] * B[(l_n*24)+12];
    C[(l_n*56)+24] += A[132] * B[(l_n*24)+12];
    C[(l_n*56)+41] += A[133] * B[(l_n*24)+12];
    C[(l_n*56)+42] += A[134] * B[(l_n*24)+12];
    C[(l_n*56)+43] += A[135] * B[(l_n*24)+12];
    C[(l_n*56)+44] += A[136] * B[(l_n*24)+12];
    C[(l_n*56)+45] += A[137] * B[(l_n*24)+12];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b13 = _mm256_broadcast_sd(&B[(l_n*24)+13]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b13 = _mm_loaddup_pd(&B[(l_n*24)+13]);
#endif
    __m128d c13_0 = _mm_load_sd(&C[(l_n*56)+20]);
    __m128d a13_0 = _mm_load_sd(&A[138]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, b13));
#endif
    _mm_store_sd(&C[(l_n*56)+20], c13_0);
    __m128d c13_1 = _mm_loadu_pd(&C[(l_n*56)+22]);
    __m128d a13_1 = _mm_loadu_pd(&A[139]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_1 = _mm_add_pd(c13_1, _mm_mul_pd(a13_1, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_1 = _mm_add_pd(c13_1, _mm_mul_pd(a13_1, b13));
#endif
    _mm_storeu_pd(&C[(l_n*56)+22], c13_1);
    __m128d c13_3 = _mm_load_sd(&C[(l_n*56)+24]);
    __m128d a13_3 = _mm_load_sd(&A[141]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_3 = _mm_add_sd(c13_3, _mm_mul_sd(a13_3, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_3 = _mm_add_sd(c13_3, _mm_mul_sd(a13_3, b13));
#endif
    _mm_store_sd(&C[(l_n*56)+24], c13_3);
    __m128d c13_4 = _mm_load_sd(&C[(l_n*56)+41]);
    __m128d a13_4 = _mm_load_sd(&A[142]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_4 = _mm_add_sd(c13_4, _mm_mul_sd(a13_4, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_4 = _mm_add_sd(c13_4, _mm_mul_sd(a13_4, b13));
#endif
    _mm_store_sd(&C[(l_n*56)+41], c13_4);
    __m128d c13_5 = _mm_loadu_pd(&C[(l_n*56)+43]);
    __m128d a13_5 = _mm_loadu_pd(&A[143]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_5 = _mm_add_pd(c13_5, _mm_mul_pd(a13_5, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_5 = _mm_add_pd(c13_5, _mm_mul_pd(a13_5, b13));
#endif
    _mm_storeu_pd(&C[(l_n*56)+43], c13_5);
    __m128d c13_7 = _mm_load_sd(&C[(l_n*56)+45]);
    __m128d a13_7 = _mm_load_sd(&A[145]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_7 = _mm_add_sd(c13_7, _mm_mul_sd(a13_7, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_7 = _mm_add_sd(c13_7, _mm_mul_sd(a13_7, b13));
#endif
    _mm_store_sd(&C[(l_n*56)+45], c13_7);
#else
    C[(l_n*56)+20] += A[138] * B[(l_n*24)+13];
    C[(l_n*56)+22] += A[139] * B[(l_n*24)+13];
    C[(l_n*56)+23] += A[140] * B[(l_n*24)+13];
    C[(l_n*56)+24] += A[141] * B[(l_n*24)+13];
    C[(l_n*56)+41] += A[142] * B[(l_n*24)+13];
    C[(l_n*56)+43] += A[143] * B[(l_n*24)+13];
    C[(l_n*56)+44] += A[144] * B[(l_n*24)+13];
    C[(l_n*56)+45] += A[145] * B[(l_n*24)+13];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b14 = _mm256_broadcast_sd(&B[(l_n*24)+14]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b14 = _mm_loaddup_pd(&B[(l_n*24)+14]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c14_0 = _mm256_loadu_pd(&C[(l_n*56)+20]);
    __m256d a14_0 = _mm256_loadu_pd(&A[146]);
    c14_0 = _mm256_add_pd(c14_0, _mm256_mul_pd(a14_0, b14));
    _mm256_storeu_pd(&C[(l_n*56)+20], c14_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c14_0 = _mm_loadu_pd(&C[(l_n*56)+20]);
    __m128d a14_0 = _mm_loadu_pd(&A[146]);
    c14_0 = _mm_add_pd(c14_0, _mm_mul_pd(a14_0, b14));
    _mm_storeu_pd(&C[(l_n*56)+20], c14_0);
    __m128d c14_2 = _mm_loadu_pd(&C[(l_n*56)+22]);
    __m128d a14_2 = _mm_loadu_pd(&A[148]);
    c14_2 = _mm_add_pd(c14_2, _mm_mul_pd(a14_2, b14));
    _mm_storeu_pd(&C[(l_n*56)+22], c14_2);
#endif
    __m128d c14_4 = _mm_load_sd(&C[(l_n*56)+24]);
    __m128d a14_4 = _mm_load_sd(&A[150]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_4 = _mm_add_sd(c14_4, _mm_mul_sd(a14_4, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_4 = _mm_add_sd(c14_4, _mm_mul_sd(a14_4, b14));
#endif
    _mm_store_sd(&C[(l_n*56)+24], c14_4);
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c14_5 = _mm256_loadu_pd(&C[(l_n*56)+41]);
    __m256d a14_5 = _mm256_loadu_pd(&A[151]);
    c14_5 = _mm256_add_pd(c14_5, _mm256_mul_pd(a14_5, b14));
    _mm256_storeu_pd(&C[(l_n*56)+41], c14_5);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c14_5 = _mm_loadu_pd(&C[(l_n*56)+41]);
    __m128d a14_5 = _mm_loadu_pd(&A[151]);
    c14_5 = _mm_add_pd(c14_5, _mm_mul_pd(a14_5, b14));
    _mm_storeu_pd(&C[(l_n*56)+41], c14_5);
    __m128d c14_7 = _mm_loadu_pd(&C[(l_n*56)+43]);
    __m128d a14_7 = _mm_loadu_pd(&A[153]);
    c14_7 = _mm_add_pd(c14_7, _mm_mul_pd(a14_7, b14));
    _mm_storeu_pd(&C[(l_n*56)+43], c14_7);
#endif
    __m128d c14_9 = _mm_load_sd(&C[(l_n*56)+45]);
    __m128d a14_9 = _mm_load_sd(&A[155]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_9 = _mm_add_sd(c14_9, _mm_mul_sd(a14_9, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_9 = _mm_add_sd(c14_9, _mm_mul_sd(a14_9, b14));
#endif
    _mm_store_sd(&C[(l_n*56)+45], c14_9);
#else
    C[(l_n*56)+20] += A[146] * B[(l_n*24)+14];
    C[(l_n*56)+21] += A[147] * B[(l_n*24)+14];
    C[(l_n*56)+22] += A[148] * B[(l_n*24)+14];
    C[(l_n*56)+23] += A[149] * B[(l_n*24)+14];
    C[(l_n*56)+24] += A[150] * B[(l_n*24)+14];
    C[(l_n*56)+41] += A[151] * B[(l_n*24)+14];
    C[(l_n*56)+42] += A[152] * B[(l_n*24)+14];
    C[(l_n*56)+43] += A[153] * B[(l_n*24)+14];
    C[(l_n*56)+44] += A[154] * B[(l_n*24)+14];
    C[(l_n*56)+45] += A[155] * B[(l_n*24)+14];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b15 = _mm256_broadcast_sd(&B[(l_n*24)+15]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b15 = _mm_loaddup_pd(&B[(l_n*24)+15]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c15_0 = _mm256_loadu_pd(&C[(l_n*56)+35]);
    __m256d a15_0 = _mm256_loadu_pd(&A[156]);
    c15_0 = _mm256_add_pd(c15_0, _mm256_mul_pd(a15_0, b15));
    _mm256_storeu_pd(&C[(l_n*56)+35], c15_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c15_0 = _mm_loadu_pd(&C[(l_n*56)+35]);
    __m128d a15_0 = _mm_loadu_pd(&A[156]);
    c15_0 = _mm_add_pd(c15_0, _mm_mul_pd(a15_0, b15));
    _mm_storeu_pd(&C[(l_n*56)+35], c15_0);
    __m128d c15_2 = _mm_loadu_pd(&C[(l_n*56)+37]);
    __m128d a15_2 = _mm_loadu_pd(&A[158]);
    c15_2 = _mm_add_pd(c15_2, _mm_mul_pd(a15_2, b15));
    _mm_storeu_pd(&C[(l_n*56)+37], c15_2);
#endif
    __m128d c15_4 = _mm_loadu_pd(&C[(l_n*56)+39]);
    __m128d a15_4 = _mm_loadu_pd(&A[160]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_4 = _mm_add_pd(c15_4, _mm_mul_pd(a15_4, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_4 = _mm_add_pd(c15_4, _mm_mul_pd(a15_4, b15));
#endif
    _mm_storeu_pd(&C[(l_n*56)+39], c15_4);
#else
    C[(l_n*56)+35] += A[156] * B[(l_n*24)+15];
    C[(l_n*56)+36] += A[157] * B[(l_n*24)+15];
    C[(l_n*56)+37] += A[158] * B[(l_n*24)+15];
    C[(l_n*56)+38] += A[159] * B[(l_n*24)+15];
    C[(l_n*56)+39] += A[160] * B[(l_n*24)+15];
    C[(l_n*56)+40] += A[161] * B[(l_n*24)+15];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b16 = _mm256_broadcast_sd(&B[(l_n*24)+16]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b16 = _mm_loaddup_pd(&B[(l_n*24)+16]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c16_0 = _mm256_loadu_pd(&C[(l_n*56)+35]);
    __m256d a16_0 = _mm256_loadu_pd(&A[162]);
    c16_0 = _mm256_add_pd(c16_0, _mm256_mul_pd(a16_0, b16));
    _mm256_storeu_pd(&C[(l_n*56)+35], c16_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c16_0 = _mm_loadu_pd(&C[(l_n*56)+35]);
    __m128d a16_0 = _mm_loadu_pd(&A[162]);
    c16_0 = _mm_add_pd(c16_0, _mm_mul_pd(a16_0, b16));
    _mm_storeu_pd(&C[(l_n*56)+35], c16_0);
    __m128d c16_2 = _mm_loadu_pd(&C[(l_n*56)+37]);
    __m128d a16_2 = _mm_loadu_pd(&A[164]);
    c16_2 = _mm_add_pd(c16_2, _mm_mul_pd(a16_2, b16));
    _mm_storeu_pd(&C[(l_n*56)+37], c16_2);
#endif
    __m128d c16_4 = _mm_loadu_pd(&C[(l_n*56)+39]);
    __m128d a16_4 = _mm_loadu_pd(&A[166]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_4 = _mm_add_pd(c16_4, _mm_mul_pd(a16_4, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_4 = _mm_add_pd(c16_4, _mm_mul_pd(a16_4, b16));
#endif
    _mm_storeu_pd(&C[(l_n*56)+39], c16_4);
#else
    C[(l_n*56)+35] += A[162] * B[(l_n*24)+16];
    C[(l_n*56)+36] += A[163] * B[(l_n*24)+16];
    C[(l_n*56)+37] += A[164] * B[(l_n*24)+16];
    C[(l_n*56)+38] += A[165] * B[(l_n*24)+16];
    C[(l_n*56)+39] += A[166] * B[(l_n*24)+16];
    C[(l_n*56)+40] += A[167] * B[(l_n*24)+16];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b17 = _mm256_broadcast_sd(&B[(l_n*24)+17]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b17 = _mm_loaddup_pd(&B[(l_n*24)+17]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c17_0 = _mm256_loadu_pd(&C[(l_n*56)+35]);
    __m256d a17_0 = _mm256_loadu_pd(&A[168]);
    c17_0 = _mm256_add_pd(c17_0, _mm256_mul_pd(a17_0, b17));
    _mm256_storeu_pd(&C[(l_n*56)+35], c17_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c17_0 = _mm_loadu_pd(&C[(l_n*56)+35]);
    __m128d a17_0 = _mm_loadu_pd(&A[168]);
    c17_0 = _mm_add_pd(c17_0, _mm_mul_pd(a17_0, b17));
    _mm_storeu_pd(&C[(l_n*56)+35], c17_0);
    __m128d c17_2 = _mm_loadu_pd(&C[(l_n*56)+37]);
    __m128d a17_2 = _mm_loadu_pd(&A[170]);
    c17_2 = _mm_add_pd(c17_2, _mm_mul_pd(a17_2, b17));
    _mm_storeu_pd(&C[(l_n*56)+37], c17_2);
#endif
    __m128d c17_4 = _mm_loadu_pd(&C[(l_n*56)+39]);
    __m128d a17_4 = _mm_loadu_pd(&A[172]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_4 = _mm_add_pd(c17_4, _mm_mul_pd(a17_4, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_4 = _mm_add_pd(c17_4, _mm_mul_pd(a17_4, b17));
#endif
    _mm_storeu_pd(&C[(l_n*56)+39], c17_4);
#else
    C[(l_n*56)+35] += A[168] * B[(l_n*24)+17];
    C[(l_n*56)+36] += A[169] * B[(l_n*24)+17];
    C[(l_n*56)+37] += A[170] * B[(l_n*24)+17];
    C[(l_n*56)+38] += A[171] * B[(l_n*24)+17];
    C[(l_n*56)+39] += A[172] * B[(l_n*24)+17];
    C[(l_n*56)+40] += A[173] * B[(l_n*24)+17];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b18 = _mm256_broadcast_sd(&B[(l_n*24)+18]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b18 = _mm_loaddup_pd(&B[(l_n*24)+18]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c18_0 = _mm256_loadu_pd(&C[(l_n*56)+35]);
    __m256d a18_0 = _mm256_loadu_pd(&A[174]);
    c18_0 = _mm256_add_pd(c18_0, _mm256_mul_pd(a18_0, b18));
    _mm256_storeu_pd(&C[(l_n*56)+35], c18_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c18_0 = _mm_loadu_pd(&C[(l_n*56)+35]);
    __m128d a18_0 = _mm_loadu_pd(&A[174]);
    c18_0 = _mm_add_pd(c18_0, _mm_mul_pd(a18_0, b18));
    _mm_storeu_pd(&C[(l_n*56)+35], c18_0);
    __m128d c18_2 = _mm_loadu_pd(&C[(l_n*56)+37]);
    __m128d a18_2 = _mm_loadu_pd(&A[176]);
    c18_2 = _mm_add_pd(c18_2, _mm_mul_pd(a18_2, b18));
    _mm_storeu_pd(&C[(l_n*56)+37], c18_2);
#endif
    __m128d c18_4 = _mm_loadu_pd(&C[(l_n*56)+39]);
    __m128d a18_4 = _mm_loadu_pd(&A[178]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_4 = _mm_add_pd(c18_4, _mm_mul_pd(a18_4, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_4 = _mm_add_pd(c18_4, _mm_mul_pd(a18_4, b18));
#endif
    _mm_storeu_pd(&C[(l_n*56)+39], c18_4);
#else
    C[(l_n*56)+35] += A[174] * B[(l_n*24)+18];
    C[(l_n*56)+36] += A[175] * B[(l_n*24)+18];
    C[(l_n*56)+37] += A[176] * B[(l_n*24)+18];
    C[(l_n*56)+38] += A[177] * B[(l_n*24)+18];
    C[(l_n*56)+39] += A[178] * B[(l_n*24)+18];
    C[(l_n*56)+40] += A[179] * B[(l_n*24)+18];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b19 = _mm256_broadcast_sd(&B[(l_n*24)+19]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b19 = _mm_loaddup_pd(&B[(l_n*24)+19]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c19_0 = _mm256_loadu_pd(&C[(l_n*56)+35]);
    __m256d a19_0 = _mm256_loadu_pd(&A[180]);
    c19_0 = _mm256_add_pd(c19_0, _mm256_mul_pd(a19_0, b19));
    _mm256_storeu_pd(&C[(l_n*56)+35], c19_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c19_0 = _mm_loadu_pd(&C[(l_n*56)+35]);
    __m128d a19_0 = _mm_loadu_pd(&A[180]);
    c19_0 = _mm_add_pd(c19_0, _mm_mul_pd(a19_0, b19));
    _mm_storeu_pd(&C[(l_n*56)+35], c19_0);
    __m128d c19_2 = _mm_loadu_pd(&C[(l_n*56)+37]);
    __m128d a19_2 = _mm_loadu_pd(&A[182]);
    c19_2 = _mm_add_pd(c19_2, _mm_mul_pd(a19_2, b19));
    _mm_storeu_pd(&C[(l_n*56)+37], c19_2);
#endif
    __m128d c19_4 = _mm_loadu_pd(&C[(l_n*56)+39]);
    __m128d a19_4 = _mm_loadu_pd(&A[184]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_4 = _mm_add_pd(c19_4, _mm_mul_pd(a19_4, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_4 = _mm_add_pd(c19_4, _mm_mul_pd(a19_4, b19));
#endif
    _mm_storeu_pd(&C[(l_n*56)+39], c19_4);
#else
    C[(l_n*56)+35] += A[180] * B[(l_n*24)+19];
    C[(l_n*56)+36] += A[181] * B[(l_n*24)+19];
    C[(l_n*56)+37] += A[182] * B[(l_n*24)+19];
    C[(l_n*56)+38] += A[183] * B[(l_n*24)+19];
    C[(l_n*56)+39] += A[184] * B[(l_n*24)+19];
    C[(l_n*56)+40] += A[185] * B[(l_n*24)+19];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b20 = _mm256_broadcast_sd(&B[(l_n*24)+20]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b20 = _mm_loaddup_pd(&B[(l_n*24)+20]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c20_0 = _mm256_loadu_pd(&C[(l_n*56)+35]);
    __m256d a20_0 = _mm256_loadu_pd(&A[186]);
    c20_0 = _mm256_add_pd(c20_0, _mm256_mul_pd(a20_0, b20));
    _mm256_storeu_pd(&C[(l_n*56)+35], c20_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c20_0 = _mm_loadu_pd(&C[(l_n*56)+35]);
    __m128d a20_0 = _mm_loadu_pd(&A[186]);
    c20_0 = _mm_add_pd(c20_0, _mm_mul_pd(a20_0, b20));
    _mm_storeu_pd(&C[(l_n*56)+35], c20_0);
    __m128d c20_2 = _mm_loadu_pd(&C[(l_n*56)+37]);
    __m128d a20_2 = _mm_loadu_pd(&A[188]);
    c20_2 = _mm_add_pd(c20_2, _mm_mul_pd(a20_2, b20));
    _mm_storeu_pd(&C[(l_n*56)+37], c20_2);
#endif
    __m128d c20_4 = _mm_loadu_pd(&C[(l_n*56)+39]);
    __m128d a20_4 = _mm_loadu_pd(&A[190]);
#if defined(__SSE3__) && defined(__AVX__)
    c20_4 = _mm_add_pd(c20_4, _mm_mul_pd(a20_4, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c20_4 = _mm_add_pd(c20_4, _mm_mul_pd(a20_4, b20));
#endif
    _mm_storeu_pd(&C[(l_n*56)+39], c20_4);
#else
    C[(l_n*56)+35] += A[186] * B[(l_n*24)+20];
    C[(l_n*56)+36] += A[187] * B[(l_n*24)+20];
    C[(l_n*56)+37] += A[188] * B[(l_n*24)+20];
    C[(l_n*56)+38] += A[189] * B[(l_n*24)+20];
    C[(l_n*56)+39] += A[190] * B[(l_n*24)+20];
    C[(l_n*56)+40] += A[191] * B[(l_n*24)+20];
#endif

  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 3456;
#endif
}

void libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vmovapd 0(%%rdx), %%ymm4\n\t"
                       "vmovapd 32(%%rdx), %%ymm5\n\t"
                       "vmovapd 64(%%rdx), %%ymm6\n\t"
                       "vmovapd 96(%%rdx), %%ymm7\n\t"
                       "vmovapd 448(%%rdx), %%ymm8\n\t"
                       "vmovapd 480(%%rdx), %%ymm9\n\t"
                       "vmovapd 512(%%rdx), %%ymm10\n\t"
                       "vmovapd 544(%%rdx), %%ymm11\n\t"
                       "vmovapd 896(%%rdx), %%ymm12\n\t"
                       "vmovapd 928(%%rdx), %%ymm13\n\t"
                       "vmovapd 960(%%rdx), %%ymm14\n\t"
                       "vmovapd 992(%%rdx), %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 384(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 392(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 400(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 216(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 408(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 224(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 416(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 232(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 424(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 240(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 432(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 248(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 440(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 256(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 264(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 456(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 272(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 464(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 280(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 472(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 480(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 296(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 488(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 304(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 496(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 312(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 504(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 320(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 512(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 328(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 520(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 336(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 528(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 344(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 536(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 352(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 544(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm4, 0(%%rdx)\n\t"
                       "vmovapd %%ymm5, 32(%%rdx)\n\t"
                       "vmovapd %%ymm6, 64(%%rdx)\n\t"
                       "vmovapd %%ymm7, 96(%%rdx)\n\t"
                       "vmovapd %%ymm8, 448(%%rdx)\n\t"
                       "vmovapd %%ymm9, 480(%%rdx)\n\t"
                       "vmovapd %%ymm10, 512(%%rdx)\n\t"
                       "vmovapd %%ymm11, 544(%%rdx)\n\t"
                       "vmovapd %%ymm12, 896(%%rdx)\n\t"
                       "vmovapd %%ymm13, 928(%%rdx)\n\t"
                       "vmovapd %%ymm14, 960(%%rdx)\n\t"
                       "vmovapd %%ymm15, 992(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $9280, %%rdi\n\t"
                       "cmpq $32, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $32, %%r10\n\t"
                       "34:\n\t"
                       "addq $12, %%r10\n\t"
                       "vmovapd 0(%%rdx), %%ymm7\n\t"
                       "vmovapd 32(%%rdx), %%ymm8\n\t"
                       "vmovapd 64(%%rdx), %%ymm9\n\t"
                       "vmovapd 448(%%rdx), %%ymm10\n\t"
                       "vmovapd 480(%%rdx), %%ymm11\n\t"
                       "vmovapd 512(%%rdx), %%ymm12\n\t"
                       "vmovapd 896(%%rdx), %%ymm13\n\t"
                       "vmovapd 928(%%rdx), %%ymm14\n\t"
                       "vmovapd 960(%%rdx), %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 384(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 392(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 400(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 216(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 408(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 224(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 416(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 232(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 424(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 240(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 432(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 248(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 440(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 256(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 264(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 456(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 272(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 464(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 280(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 472(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 480(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 296(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 488(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 304(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 496(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 312(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 504(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 320(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 512(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 328(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 520(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 336(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 528(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 344(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 536(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 352(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 544(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm7, 0(%%rdx)\n\t"
                       "vmovapd %%ymm8, 32(%%rdx)\n\t"
                       "vmovapd %%ymm9, 64(%%rdx)\n\t"
                       "vmovapd %%ymm10, 448(%%rdx)\n\t"
                       "vmovapd %%ymm11, 480(%%rdx)\n\t"
                       "vmovapd %%ymm12, 512(%%rdx)\n\t"
                       "vmovapd %%ymm13, 896(%%rdx)\n\t"
                       "vmovapd %%ymm14, 928(%%rdx)\n\t"
                       "vmovapd %%ymm15, 960(%%rdx)\n\t"
                       "addq $96, %%rdx\n\t"
                       "subq $9312, %%rdi\n\t"
                       "cmpq $56, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $896, %%rdx\n\t"
                       "addq $576, %%rsi\n\t"
                       "subq $448, %%rdi\n\t"
                       "cmpq $9, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 21168;
#endif
}

void libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vxorpd %%ymm4, %%ymm4, %%ymm4\n\t"
                       "vxorpd %%ymm5, %%ymm5, %%ymm5\n\t"
                       "vxorpd %%ymm6, %%ymm6, %%ymm6\n\t"
                       "vxorpd %%ymm7, %%ymm7, %%ymm7\n\t"
                       "vxorpd %%ymm8, %%ymm8, %%ymm8\n\t"
                       "vxorpd %%ymm9, %%ymm9, %%ymm9\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 384(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 392(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 400(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 216(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 408(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 224(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 416(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 232(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 424(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 240(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 432(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 248(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 440(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 256(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 264(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 456(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 272(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 464(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 280(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 472(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 480(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 296(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 488(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 304(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 496(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 312(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 504(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 320(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 512(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 328(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 520(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 336(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 528(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 344(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 536(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 352(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 544(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm4, 0(%%rdx)\n\t"
                       "vmovapd %%ymm5, 32(%%rdx)\n\t"
                       "vmovapd %%ymm6, 64(%%rdx)\n\t"
                       "vmovapd %%ymm7, 96(%%rdx)\n\t"
                       "vmovapd %%ymm8, 192(%%rdx)\n\t"
                       "vmovapd %%ymm9, 224(%%rdx)\n\t"
                       "vmovapd %%ymm10, 256(%%rdx)\n\t"
                       "vmovapd %%ymm11, 288(%%rdx)\n\t"
                       "vmovapd %%ymm12, 384(%%rdx)\n\t"
                       "vmovapd %%ymm13, 416(%%rdx)\n\t"
                       "vmovapd %%ymm14, 448(%%rdx)\n\t"
                       "vmovapd %%ymm15, 480(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $3904, %%rdi\n\t"
                       "cmpq $16, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $16, %%r10\n\t"
                       "34:\n\t"
                       "addq $8, %%r10\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 384(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 392(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 400(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 216(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 408(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 224(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 416(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 232(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 424(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 240(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 432(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 248(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 440(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 256(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 264(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 456(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 272(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 464(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 280(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 472(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 480(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 296(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 488(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 304(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 496(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 312(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 504(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 320(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 512(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 328(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 520(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 336(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 528(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 344(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 536(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 352(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 544(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm10, 0(%%rdx)\n\t"
                       "vmovapd %%ymm11, 32(%%rdx)\n\t"
                       "vmovapd %%ymm12, 192(%%rdx)\n\t"
                       "vmovapd %%ymm13, 224(%%rdx)\n\t"
                       "vmovapd %%ymm14, 384(%%rdx)\n\t"
                       "vmovapd %%ymm15, 416(%%rdx)\n\t"
                       "addq $64, %%rdx\n\t"
                       "subq $3968, %%rdi\n\t"
                       "cmpq $24, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $384, %%rdx\n\t"
                       "addq $576, %%rsi\n\t"
                       "subq $192, %%rdi\n\t"
                       "cmpq $9, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 9072;
#endif
}

void libxsmmsparse_b8016671674e78f7e2ec1c02d0b31df8_m21_n9_k21_ldA0_ldB24_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_n = 0;
  #pragma nounroll_and_jam
  for ( l_n = 0; l_n < 9; l_n++) {
    unsigned int l_m = 0;
   #pragma simd
    for ( l_m = 0; l_m < 21; l_m++) {
      C[(l_n*21)+l_m] = 0.0;
    }
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b0 = _mm256_broadcast_sd(&B[(l_n*24)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b0 = _mm_loaddup_pd(&B[(l_n*24)+0]);
#endif
    __m128d c0_0 = _mm_load_sd(&C[(l_n*21)+0]);
    __m128d a0_0 = _mm_load_sd(&A[0]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, b0));
#endif
    _mm_store_sd(&C[(l_n*21)+0], c0_0);
#else
    C[(l_n*21)+0] += A[0] * B[(l_n*24)+0];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b1 = _mm256_broadcast_sd(&B[(l_n*24)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b1 = _mm_loaddup_pd(&B[(l_n*24)+1]);
#endif
    __m128d c1_0 = _mm_load_sd(&C[(l_n*21)+1]);
    __m128d a1_0 = _mm_load_sd(&A[1]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_0 = _mm_add_sd(c1_0, _mm_mul_sd(a1_0, b1));
#endif
    _mm_store_sd(&C[(l_n*21)+1], c1_0);
#else
    C[(l_n*21)+1] += A[1] * B[(l_n*24)+1];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b2 = _mm256_broadcast_sd(&B[(l_n*24)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b2 = _mm_loaddup_pd(&B[(l_n*24)+2]);
#endif
    __m128d c2_0 = _mm_load_sd(&C[(l_n*21)+2]);
    __m128d a2_0 = _mm_load_sd(&A[2]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_0 = _mm_add_sd(c2_0, _mm_mul_sd(a2_0, b2));
#endif
    _mm_store_sd(&C[(l_n*21)+2], c2_0);
#else
    C[(l_n*21)+2] += A[2] * B[(l_n*24)+2];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b3 = _mm256_broadcast_sd(&B[(l_n*24)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b3 = _mm_loaddup_pd(&B[(l_n*24)+3]);
#endif
    __m128d c3_0 = _mm_load_sd(&C[(l_n*21)+3]);
    __m128d a3_0 = _mm_load_sd(&A[3]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_0 = _mm_add_sd(c3_0, _mm_mul_sd(a3_0, b3));
#endif
    _mm_store_sd(&C[(l_n*21)+3], c3_0);
#else
    C[(l_n*21)+3] += A[3] * B[(l_n*24)+3];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b4 = _mm256_broadcast_sd(&B[(l_n*24)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b4 = _mm_loaddup_pd(&B[(l_n*24)+4]);
#endif
    __m128d c4_0 = _mm_load_sd(&C[(l_n*21)+4]);
    __m128d a4_0 = _mm_load_sd(&A[4]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_0 = _mm_add_sd(c4_0, _mm_mul_sd(a4_0, b4));
#endif
    _mm_store_sd(&C[(l_n*21)+4], c4_0);
#else
    C[(l_n*21)+4] += A[4] * B[(l_n*24)+4];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b5 = _mm256_broadcast_sd(&B[(l_n*24)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b5 = _mm_loaddup_pd(&B[(l_n*24)+5]);
#endif
    __m128d c5_0 = _mm_load_sd(&C[(l_n*21)+5]);
    __m128d a5_0 = _mm_load_sd(&A[5]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_0 = _mm_add_sd(c5_0, _mm_mul_sd(a5_0, b5));
#endif
    _mm_store_sd(&C[(l_n*21)+5], c5_0);
#else
    C[(l_n*21)+5] += A[5] * B[(l_n*24)+5];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b6 = _mm256_broadcast_sd(&B[(l_n*24)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b6 = _mm_loaddup_pd(&B[(l_n*24)+6]);
#endif
    __m128d c6_0 = _mm_load_sd(&C[(l_n*21)+6]);
    __m128d a6_0 = _mm_load_sd(&A[6]);
#if defined(__SSE3__) && defined(__AVX__)
    c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, _mm256_castpd256_pd128(b6)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c6_0 = _mm_add_sd(c6_0, _mm_mul_sd(a6_0, b6));
#endif
    _mm_store_sd(&C[(l_n*21)+6], c6_0);
#else
    C[(l_n*21)+6] += A[6] * B[(l_n*24)+6];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b7 = _mm256_broadcast_sd(&B[(l_n*24)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b7 = _mm_loaddup_pd(&B[(l_n*24)+7]);
#endif
    __m128d c7_0 = _mm_load_sd(&C[(l_n*21)+7]);
    __m128d a7_0 = _mm_load_sd(&A[7]);
#if defined(__SSE3__) && defined(__AVX__)
    c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, _mm256_castpd256_pd128(b7)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c7_0 = _mm_add_sd(c7_0, _mm_mul_sd(a7_0, b7));
#endif
    _mm_store_sd(&C[(l_n*21)+7], c7_0);
#else
    C[(l_n*21)+7] += A[7] * B[(l_n*24)+7];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b8 = _mm256_broadcast_sd(&B[(l_n*24)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b8 = _mm_loaddup_pd(&B[(l_n*24)+8]);
#endif
    __m128d c8_0 = _mm_load_sd(&C[(l_n*21)+8]);
    __m128d a8_0 = _mm_load_sd(&A[8]);
#if defined(__SSE3__) && defined(__AVX__)
    c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, _mm256_castpd256_pd128(b8)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c8_0 = _mm_add_sd(c8_0, _mm_mul_sd(a8_0, b8));
#endif
    _mm_store_sd(&C[(l_n*21)+8], c8_0);
#else
    C[(l_n*21)+8] += A[8] * B[(l_n*24)+8];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b9 = _mm256_broadcast_sd(&B[(l_n*24)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b9 = _mm_loaddup_pd(&B[(l_n*24)+9]);
#endif
    __m128d c9_0 = _mm_load_sd(&C[(l_n*21)+9]);
    __m128d a9_0 = _mm_load_sd(&A[9]);
#if defined(__SSE3__) && defined(__AVX__)
    c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, _mm256_castpd256_pd128(b9)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c9_0 = _mm_add_sd(c9_0, _mm_mul_sd(a9_0, b9));
#endif
    _mm_store_sd(&C[(l_n*21)+9], c9_0);
#else
    C[(l_n*21)+9] += A[9] * B[(l_n*24)+9];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b10 = _mm256_broadcast_sd(&B[(l_n*24)+10]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b10 = _mm_loaddup_pd(&B[(l_n*24)+10]);
#endif
    __m128d c10_0 = _mm_load_sd(&C[(l_n*21)+10]);
    __m128d a10_0 = _mm_load_sd(&A[10]);
#if defined(__SSE3__) && defined(__AVX__)
    c10_0 = _mm_add_sd(c10_0, _mm_mul_sd(a10_0, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c10_0 = _mm_add_sd(c10_0, _mm_mul_sd(a10_0, b10));
#endif
    _mm_store_sd(&C[(l_n*21)+10], c10_0);
#else
    C[(l_n*21)+10] += A[10] * B[(l_n*24)+10];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b11 = _mm256_broadcast_sd(&B[(l_n*24)+11]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b11 = _mm_loaddup_pd(&B[(l_n*24)+11]);
#endif
    __m128d c11_0 = _mm_load_sd(&C[(l_n*21)+11]);
    __m128d a11_0 = _mm_load_sd(&A[11]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_0 = _mm_add_sd(c11_0, _mm_mul_sd(a11_0, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_0 = _mm_add_sd(c11_0, _mm_mul_sd(a11_0, b11));
#endif
    _mm_store_sd(&C[(l_n*21)+11], c11_0);
#else
    C[(l_n*21)+11] += A[11] * B[(l_n*24)+11];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b12 = _mm256_broadcast_sd(&B[(l_n*24)+12]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b12 = _mm_loaddup_pd(&B[(l_n*24)+12]);
#endif
    __m128d c12_0 = _mm_load_sd(&C[(l_n*21)+12]);
    __m128d a12_0 = _mm_load_sd(&A[12]);
#if defined(__SSE3__) && defined(__AVX__)
    c12_0 = _mm_add_sd(c12_0, _mm_mul_sd(a12_0, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c12_0 = _mm_add_sd(c12_0, _mm_mul_sd(a12_0, b12));
#endif
    _mm_store_sd(&C[(l_n*21)+12], c12_0);
#else
    C[(l_n*21)+12] += A[12] * B[(l_n*24)+12];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b13 = _mm256_broadcast_sd(&B[(l_n*24)+13]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b13 = _mm_loaddup_pd(&B[(l_n*24)+13]);
#endif
    __m128d c13_0 = _mm_load_sd(&C[(l_n*21)+13]);
    __m128d a13_0 = _mm_load_sd(&A[13]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, b13));
#endif
    _mm_store_sd(&C[(l_n*21)+13], c13_0);
#else
    C[(l_n*21)+13] += A[13] * B[(l_n*24)+13];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b14 = _mm256_broadcast_sd(&B[(l_n*24)+14]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b14 = _mm_loaddup_pd(&B[(l_n*24)+14]);
#endif
    __m128d c14_0 = _mm_load_sd(&C[(l_n*21)+14]);
    __m128d a14_0 = _mm_load_sd(&A[14]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_0 = _mm_add_sd(c14_0, _mm_mul_sd(a14_0, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_0 = _mm_add_sd(c14_0, _mm_mul_sd(a14_0, b14));
#endif
    _mm_store_sd(&C[(l_n*21)+14], c14_0);
#else
    C[(l_n*21)+14] += A[14] * B[(l_n*24)+14];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b15 = _mm256_broadcast_sd(&B[(l_n*24)+15]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b15 = _mm_loaddup_pd(&B[(l_n*24)+15]);
#endif
    __m128d c15_0 = _mm_load_sd(&C[(l_n*21)+15]);
    __m128d a15_0 = _mm_load_sd(&A[15]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_0 = _mm_add_sd(c15_0, _mm_mul_sd(a15_0, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_0 = _mm_add_sd(c15_0, _mm_mul_sd(a15_0, b15));
#endif
    _mm_store_sd(&C[(l_n*21)+15], c15_0);
#else
    C[(l_n*21)+15] += A[15] * B[(l_n*24)+15];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b16 = _mm256_broadcast_sd(&B[(l_n*24)+16]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b16 = _mm_loaddup_pd(&B[(l_n*24)+16]);
#endif
    __m128d c16_0 = _mm_load_sd(&C[(l_n*21)+16]);
    __m128d a16_0 = _mm_load_sd(&A[16]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_0 = _mm_add_sd(c16_0, _mm_mul_sd(a16_0, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_0 = _mm_add_sd(c16_0, _mm_mul_sd(a16_0, b16));
#endif
    _mm_store_sd(&C[(l_n*21)+16], c16_0);
#else
    C[(l_n*21)+16] += A[16] * B[(l_n*24)+16];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b17 = _mm256_broadcast_sd(&B[(l_n*24)+17]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b17 = _mm_loaddup_pd(&B[(l_n*24)+17]);
#endif
    __m128d c17_0 = _mm_load_sd(&C[(l_n*21)+17]);
    __m128d a17_0 = _mm_load_sd(&A[17]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_0 = _mm_add_sd(c17_0, _mm_mul_sd(a17_0, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_0 = _mm_add_sd(c17_0, _mm_mul_sd(a17_0, b17));
#endif
    _mm_store_sd(&C[(l_n*21)+17], c17_0);
#else
    C[(l_n*21)+17] += A[17] * B[(l_n*24)+17];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b18 = _mm256_broadcast_sd(&B[(l_n*24)+18]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b18 = _mm_loaddup_pd(&B[(l_n*24)+18]);
#endif
    __m128d c18_0 = _mm_load_sd(&C[(l_n*21)+18]);
    __m128d a18_0 = _mm_load_sd(&A[18]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_0 = _mm_add_sd(c18_0, _mm_mul_sd(a18_0, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_0 = _mm_add_sd(c18_0, _mm_mul_sd(a18_0, b18));
#endif
    _mm_store_sd(&C[(l_n*21)+18], c18_0);
#else
    C[(l_n*21)+18] += A[18] * B[(l_n*24)+18];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b19 = _mm256_broadcast_sd(&B[(l_n*24)+19]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b19 = _mm_loaddup_pd(&B[(l_n*24)+19]);
#endif
    __m128d c19_0 = _mm_load_sd(&C[(l_n*21)+19]);
    __m128d a19_0 = _mm_load_sd(&A[19]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_0 = _mm_add_sd(c19_0, _mm_mul_sd(a19_0, b19));
#endif
    _mm_store_sd(&C[(l_n*21)+19], c19_0);
#else
    C[(l_n*21)+19] += A[19] * B[(l_n*24)+19];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b20 = _mm256_broadcast_sd(&B[(l_n*24)+20]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b20 = _mm_loaddup_pd(&B[(l_n*24)+20]);
#endif
    __m128d c20_0 = _mm_load_sd(&C[(l_n*21)+20]);
    __m128d a20_0 = _mm_load_sd(&A[20]);
#if defined(__SSE3__) && defined(__AVX__)
    c20_0 = _mm_add_sd(c20_0, _mm_mul_sd(a20_0, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c20_0 = _mm_add_sd(c20_0, _mm_mul_sd(a20_0, b20));
#endif
    _mm_store_sd(&C[(l_n*21)+20], c20_0);
#else
    C[(l_n*21)+20] += A[20] * B[(l_n*24)+20];
#endif

  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 378;
#endif
}

void libxsmm_m21_n9_k9_ldA21_ldB9_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vxorpd %%ymm4, %%ymm4, %%ymm4\n\t"
                       "vxorpd %%ymm5, %%ymm5, %%ymm5\n\t"
                       "vxorpd %%ymm6, %%ymm6, %%ymm6\n\t"
                       "vxorpd %%ymm7, %%ymm7, %%ymm7\n\t"
                       "vxorpd %%ymm8, %%ymm8, %%ymm8\n\t"
                       "vxorpd %%ymm9, %%ymm9, %%ymm9\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm2\n\t"
                       "vmovupd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovupd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovupd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovupd 96(%%rdi), %%ymm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm2\n\t"
                       "vmovupd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovupd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovupd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovupd 96(%%rdi), %%ymm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm2\n\t"
                       "vmovupd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovupd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovupd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovupd 96(%%rdi), %%ymm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 168(%%rsi), %%ymm2\n\t"
                       "vmovupd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovupd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovupd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovupd 96(%%rdi), %%ymm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 176(%%rsi), %%ymm2\n\t"
                       "vmovupd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovupd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovupd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovupd 96(%%rdi), %%ymm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 184(%%rsi), %%ymm2\n\t"
                       "vmovupd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovupd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovupd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovupd 96(%%rdi), %%ymm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm2\n\t"
                       "vmovupd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovupd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovupd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovupd 96(%%rdi), %%ymm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm2\n\t"
                       "vmovupd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovupd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovupd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovupd 96(%%rdi), %%ymm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm2\n\t"
                       "vmovupd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovupd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovupd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovupd 96(%%rdi), %%ymm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovupd %%ymm4, 0(%%rdx)\n\t"
                       "vmovupd %%ymm5, 32(%%rdx)\n\t"
                       "vmovupd %%ymm6, 64(%%rdx)\n\t"
                       "vmovupd %%ymm7, 96(%%rdx)\n\t"
                       "vmovupd %%ymm8, 168(%%rdx)\n\t"
                       "vmovupd %%ymm9, 200(%%rdx)\n\t"
                       "vmovupd %%ymm10, 232(%%rdx)\n\t"
                       "vmovupd %%ymm11, 264(%%rdx)\n\t"
                       "vmovupd %%ymm12, 336(%%rdx)\n\t"
                       "vmovupd %%ymm13, 368(%%rdx)\n\t"
                       "vmovupd %%ymm14, 400(%%rdx)\n\t"
                       "vmovupd %%ymm15, 432(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $1384, %%rdi\n\t"
                       "cmpq $16, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $16, %%r10\n\t"
                       "34:\n\t"
                       "addq $4, %%r10\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 168(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 176(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 184(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovupd 0(%%rdi), %%ymm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovupd %%ymm13, 0(%%rdx)\n\t"
                       "vmovupd %%ymm14, 168(%%rdx)\n\t"
                       "vmovupd %%ymm15, 336(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $1480, %%rdi\n\t"
                       "cmpq $20, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $20, %%r10\n\t"
                       "34:\n\t"
                       "addq $1, %%r10\n\t"
                       "vxorpd %%xmm13, %%xmm13, %%xmm13\n\t"
                       "vxorpd %%xmm14, %%xmm14, %%xmm14\n\t"
                       "vxorpd %%xmm15, %%xmm15, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vmovsd 0(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm3, %%xmm0, %%xmm13\n\t"
                       "vmovsd 72(%%rsi), %%xmm1\n\t"
                       "vfmadd231sd %%xmm3, %%xmm1, %%xmm14\n\t"
                       "vmovsd 144(%%rsi), %%xmm2\n\t"
                       "vfmadd231sd %%xmm3, %%xmm2, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vmovsd 8(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm3, %%xmm0, %%xmm13\n\t"
                       "vmovsd 80(%%rsi), %%xmm1\n\t"
                       "vfmadd231sd %%xmm3, %%xmm1, %%xmm14\n\t"
                       "vmovsd 152(%%rsi), %%xmm2\n\t"
                       "vfmadd231sd %%xmm3, %%xmm2, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vmovsd 16(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm3, %%xmm0, %%xmm13\n\t"
                       "vmovsd 88(%%rsi), %%xmm1\n\t"
                       "vfmadd231sd %%xmm3, %%xmm1, %%xmm14\n\t"
                       "vmovsd 160(%%rsi), %%xmm2\n\t"
                       "vfmadd231sd %%xmm3, %%xmm2, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vmovsd 24(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm3, %%xmm0, %%xmm13\n\t"
                       "vmovsd 96(%%rsi), %%xmm1\n\t"
                       "vfmadd231sd %%xmm3, %%xmm1, %%xmm14\n\t"
                       "vmovsd 168(%%rsi), %%xmm2\n\t"
                       "vfmadd231sd %%xmm3, %%xmm2, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vmovsd 32(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm3, %%xmm0, %%xmm13\n\t"
                       "vmovsd 104(%%rsi), %%xmm1\n\t"
                       "vfmadd231sd %%xmm3, %%xmm1, %%xmm14\n\t"
                       "vmovsd 176(%%rsi), %%xmm2\n\t"
                       "vfmadd231sd %%xmm3, %%xmm2, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vmovsd 40(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm3, %%xmm0, %%xmm13\n\t"
                       "vmovsd 112(%%rsi), %%xmm1\n\t"
                       "vfmadd231sd %%xmm3, %%xmm1, %%xmm14\n\t"
                       "vmovsd 184(%%rsi), %%xmm2\n\t"
                       "vfmadd231sd %%xmm3, %%xmm2, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vmovsd 48(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm3, %%xmm0, %%xmm13\n\t"
                       "vmovsd 120(%%rsi), %%xmm1\n\t"
                       "vfmadd231sd %%xmm3, %%xmm1, %%xmm14\n\t"
                       "vmovsd 192(%%rsi), %%xmm2\n\t"
                       "vfmadd231sd %%xmm3, %%xmm2, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vmovsd 56(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm3, %%xmm0, %%xmm13\n\t"
                       "vmovsd 128(%%rsi), %%xmm1\n\t"
                       "vfmadd231sd %%xmm3, %%xmm1, %%xmm14\n\t"
                       "vmovsd 200(%%rsi), %%xmm2\n\t"
                       "vfmadd231sd %%xmm3, %%xmm2, %%xmm15\n\t"
                       "vmovsd 0(%%rdi), %%xmm3\n\t"
                       "addq $168, %%rdi\n\t"
                       "vmovsd 64(%%rsi), %%xmm0\n\t"
                       "vfmadd231sd %%xmm3, %%xmm0, %%xmm13\n\t"
                       "vmovsd 136(%%rsi), %%xmm1\n\t"
                       "vfmadd231sd %%xmm3, %%xmm1, %%xmm14\n\t"
                       "vmovsd 208(%%rsi), %%xmm2\n\t"
                       "vfmadd231sd %%xmm3, %%xmm2, %%xmm15\n\t"
                       "vmovsd %%xmm13, 0(%%rdx)\n\t"
                       "vmovsd %%xmm14, 168(%%rdx)\n\t"
                       "vmovsd %%xmm15, 336(%%rdx)\n\t"
                       "addq $8, %%rdx\n\t"
                       "subq $1504, %%rdi\n\t"
                       "cmpq $21, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $336, %%rdx\n\t"
                       "addq $216, %%rsi\n\t"
                       "subq $168, %%rdi\n\t"
                       "cmpq $9, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 3402;
#endif
}

void libxsmmsparse_38193e2bc92c2e551f139883f9580e10_m56_n9_k21_ldA0_ldB21_ldC56_alpha1_beta1_alignedA0_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_n = 0;
  #pragma nounroll_and_jam
  for ( l_n = 0; l_n < 9; l_n++) {
#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b0 = _mm256_broadcast_sd(&B[(l_n*21)+0]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b0 = _mm_loaddup_pd(&B[(l_n*21)+0]);
#endif
    __m128d c0_0 = _mm_load_sd(&C[(l_n*56)+0]);
    __m128d a0_0 = _mm_load_sd(&A[0]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_0 = _mm_add_sd(c0_0, _mm_mul_sd(a0_0, b0));
#endif
    _mm_store_sd(&C[(l_n*56)+0], c0_0);
    __m128d c0_1 = _mm_load_sd(&C[(l_n*56)+3]);
    __m128d a0_1 = _mm_load_sd(&A[1]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_1 = _mm_add_sd(c0_1, _mm_mul_sd(a0_1, b0));
#endif
    _mm_store_sd(&C[(l_n*56)+3], c0_1);
    __m128d c0_2 = _mm_load_sd(&C[(l_n*56)+9]);
    __m128d a0_2 = _mm_load_sd(&A[2]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_2 = _mm_add_sd(c0_2, _mm_mul_sd(a0_2, b0));
#endif
    _mm_store_sd(&C[(l_n*56)+9], c0_2);
    __m128d c0_3 = _mm_load_sd(&C[(l_n*56)+19]);
    __m128d a0_3 = _mm_load_sd(&A[3]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_3 = _mm_add_sd(c0_3, _mm_mul_sd(a0_3, b0));
#endif
    _mm_store_sd(&C[(l_n*56)+19], c0_3);
    __m128d c0_4 = _mm_load_sd(&C[(l_n*56)+34]);
    __m128d a0_4 = _mm_load_sd(&A[4]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_4 = _mm_add_sd(c0_4, _mm_mul_sd(a0_4, b0));
#endif
    _mm_store_sd(&C[(l_n*56)+34], c0_4);
    __m128d c0_5 = _mm_load_sd(&C[(l_n*56)+55]);
    __m128d a0_5 = _mm_load_sd(&A[5]);
#if defined(__SSE3__) && defined(__AVX__)
    c0_5 = _mm_add_sd(c0_5, _mm_mul_sd(a0_5, _mm256_castpd256_pd128(b0)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c0_5 = _mm_add_sd(c0_5, _mm_mul_sd(a0_5, b0));
#endif
    _mm_store_sd(&C[(l_n*56)+55], c0_5);
#else
    C[(l_n*56)+0] += A[0] * B[(l_n*21)+0];
    C[(l_n*56)+3] += A[1] * B[(l_n*21)+0];
    C[(l_n*56)+9] += A[2] * B[(l_n*21)+0];
    C[(l_n*56)+19] += A[3] * B[(l_n*21)+0];
    C[(l_n*56)+34] += A[4] * B[(l_n*21)+0];
    C[(l_n*56)+55] += A[5] * B[(l_n*21)+0];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b1 = _mm256_broadcast_sd(&B[(l_n*21)+1]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b1 = _mm_loaddup_pd(&B[(l_n*21)+1]);
#endif
    __m128d c1_0 = _mm_loadu_pd(&C[(l_n*56)+1]);
    __m128d a1_0 = _mm_loadu_pd(&A[6]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_0 = _mm_add_pd(c1_0, _mm_mul_pd(a1_0, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_0 = _mm_add_pd(c1_0, _mm_mul_pd(a1_0, b1));
#endif
    _mm_storeu_pd(&C[(l_n*56)+1], c1_0);
    __m128d c1_2 = _mm_loadu_pd(&C[(l_n*56)+7]);
    __m128d a1_2 = _mm_loadu_pd(&A[8]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_2 = _mm_add_pd(c1_2, _mm_mul_pd(a1_2, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_2 = _mm_add_pd(c1_2, _mm_mul_pd(a1_2, b1));
#endif
    _mm_storeu_pd(&C[(l_n*56)+7], c1_2);
    __m128d c1_4 = _mm_loadu_pd(&C[(l_n*56)+17]);
    __m128d a1_4 = _mm_loadu_pd(&A[10]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_4 = _mm_add_pd(c1_4, _mm_mul_pd(a1_4, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_4 = _mm_add_pd(c1_4, _mm_mul_pd(a1_4, b1));
#endif
    _mm_storeu_pd(&C[(l_n*56)+17], c1_4);
    __m128d c1_6 = _mm_loadu_pd(&C[(l_n*56)+32]);
    __m128d a1_6 = _mm_loadu_pd(&A[12]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_6 = _mm_add_pd(c1_6, _mm_mul_pd(a1_6, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_6 = _mm_add_pd(c1_6, _mm_mul_pd(a1_6, b1));
#endif
    _mm_storeu_pd(&C[(l_n*56)+32], c1_6);
    __m128d c1_8 = _mm_loadu_pd(&C[(l_n*56)+53]);
    __m128d a1_8 = _mm_loadu_pd(&A[14]);
#if defined(__SSE3__) && defined(__AVX__)
    c1_8 = _mm_add_pd(c1_8, _mm_mul_pd(a1_8, _mm256_castpd256_pd128(b1)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c1_8 = _mm_add_pd(c1_8, _mm_mul_pd(a1_8, b1));
#endif
    _mm_storeu_pd(&C[(l_n*56)+53], c1_8);
#else
    C[(l_n*56)+1] += A[6] * B[(l_n*21)+1];
    C[(l_n*56)+2] += A[7] * B[(l_n*21)+1];
    C[(l_n*56)+7] += A[8] * B[(l_n*21)+1];
    C[(l_n*56)+8] += A[9] * B[(l_n*21)+1];
    C[(l_n*56)+17] += A[10] * B[(l_n*21)+1];
    C[(l_n*56)+18] += A[11] * B[(l_n*21)+1];
    C[(l_n*56)+32] += A[12] * B[(l_n*21)+1];
    C[(l_n*56)+33] += A[13] * B[(l_n*21)+1];
    C[(l_n*56)+53] += A[14] * B[(l_n*21)+1];
    C[(l_n*56)+54] += A[15] * B[(l_n*21)+1];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b2 = _mm256_broadcast_sd(&B[(l_n*21)+2]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b2 = _mm_loaddup_pd(&B[(l_n*21)+2]);
#endif
    __m128d c2_0 = _mm_loadu_pd(&C[(l_n*56)+1]);
    __m128d a2_0 = _mm_loadu_pd(&A[16]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_0 = _mm_add_pd(c2_0, _mm_mul_pd(a2_0, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_0 = _mm_add_pd(c2_0, _mm_mul_pd(a2_0, b2));
#endif
    _mm_storeu_pd(&C[(l_n*56)+1], c2_0);
    __m128d c2_2 = _mm_loadu_pd(&C[(l_n*56)+7]);
    __m128d a2_2 = _mm_loadu_pd(&A[18]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_2 = _mm_add_pd(c2_2, _mm_mul_pd(a2_2, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_2 = _mm_add_pd(c2_2, _mm_mul_pd(a2_2, b2));
#endif
    _mm_storeu_pd(&C[(l_n*56)+7], c2_2);
    __m128d c2_4 = _mm_loadu_pd(&C[(l_n*56)+17]);
    __m128d a2_4 = _mm_loadu_pd(&A[20]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_4 = _mm_add_pd(c2_4, _mm_mul_pd(a2_4, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_4 = _mm_add_pd(c2_4, _mm_mul_pd(a2_4, b2));
#endif
    _mm_storeu_pd(&C[(l_n*56)+17], c2_4);
    __m128d c2_6 = _mm_loadu_pd(&C[(l_n*56)+32]);
    __m128d a2_6 = _mm_loadu_pd(&A[22]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_6 = _mm_add_pd(c2_6, _mm_mul_pd(a2_6, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_6 = _mm_add_pd(c2_6, _mm_mul_pd(a2_6, b2));
#endif
    _mm_storeu_pd(&C[(l_n*56)+32], c2_6);
    __m128d c2_8 = _mm_loadu_pd(&C[(l_n*56)+53]);
    __m128d a2_8 = _mm_loadu_pd(&A[24]);
#if defined(__SSE3__) && defined(__AVX__)
    c2_8 = _mm_add_pd(c2_8, _mm_mul_pd(a2_8, _mm256_castpd256_pd128(b2)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c2_8 = _mm_add_pd(c2_8, _mm_mul_pd(a2_8, b2));
#endif
    _mm_storeu_pd(&C[(l_n*56)+53], c2_8);
#else
    C[(l_n*56)+1] += A[16] * B[(l_n*21)+2];
    C[(l_n*56)+2] += A[17] * B[(l_n*21)+2];
    C[(l_n*56)+7] += A[18] * B[(l_n*21)+2];
    C[(l_n*56)+8] += A[19] * B[(l_n*21)+2];
    C[(l_n*56)+17] += A[20] * B[(l_n*21)+2];
    C[(l_n*56)+18] += A[21] * B[(l_n*21)+2];
    C[(l_n*56)+32] += A[22] * B[(l_n*21)+2];
    C[(l_n*56)+33] += A[23] * B[(l_n*21)+2];
    C[(l_n*56)+53] += A[24] * B[(l_n*21)+2];
    C[(l_n*56)+54] += A[25] * B[(l_n*21)+2];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b3 = _mm256_broadcast_sd(&B[(l_n*21)+3]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b3 = _mm_loaddup_pd(&B[(l_n*21)+3]);
#endif
    __m128d c3_0 = _mm_loadu_pd(&C[(l_n*56)+4]);
    __m128d a3_0 = _mm_loadu_pd(&A[26]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_0 = _mm_add_pd(c3_0, _mm_mul_pd(a3_0, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_0 = _mm_add_pd(c3_0, _mm_mul_pd(a3_0, b3));
#endif
    _mm_storeu_pd(&C[(l_n*56)+4], c3_0);
    __m128d c3_2 = _mm_load_sd(&C[(l_n*56)+6]);
    __m128d a3_2 = _mm_load_sd(&A[28]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_2 = _mm_add_sd(c3_2, _mm_mul_sd(a3_2, b3));
#endif
    _mm_store_sd(&C[(l_n*56)+6], c3_2);
    __m128d c3_3 = _mm_loadu_pd(&C[(l_n*56)+14]);
    __m128d a3_3 = _mm_loadu_pd(&A[29]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_3 = _mm_add_pd(c3_3, _mm_mul_pd(a3_3, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_3 = _mm_add_pd(c3_3, _mm_mul_pd(a3_3, b3));
#endif
    _mm_storeu_pd(&C[(l_n*56)+14], c3_3);
    __m128d c3_5 = _mm_load_sd(&C[(l_n*56)+16]);
    __m128d a3_5 = _mm_load_sd(&A[31]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_5 = _mm_add_sd(c3_5, _mm_mul_sd(a3_5, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_5 = _mm_add_sd(c3_5, _mm_mul_sd(a3_5, b3));
#endif
    _mm_store_sd(&C[(l_n*56)+16], c3_5);
    __m128d c3_6 = _mm_loadu_pd(&C[(l_n*56)+29]);
    __m128d a3_6 = _mm_loadu_pd(&A[32]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_6 = _mm_add_pd(c3_6, _mm_mul_pd(a3_6, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_6 = _mm_add_pd(c3_6, _mm_mul_pd(a3_6, b3));
#endif
    _mm_storeu_pd(&C[(l_n*56)+29], c3_6);
    __m128d c3_8 = _mm_load_sd(&C[(l_n*56)+31]);
    __m128d a3_8 = _mm_load_sd(&A[34]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_8 = _mm_add_sd(c3_8, _mm_mul_sd(a3_8, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_8 = _mm_add_sd(c3_8, _mm_mul_sd(a3_8, b3));
#endif
    _mm_store_sd(&C[(l_n*56)+31], c3_8);
    __m128d c3_9 = _mm_loadu_pd(&C[(l_n*56)+50]);
    __m128d a3_9 = _mm_loadu_pd(&A[35]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_9 = _mm_add_pd(c3_9, _mm_mul_pd(a3_9, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_9 = _mm_add_pd(c3_9, _mm_mul_pd(a3_9, b3));
#endif
    _mm_storeu_pd(&C[(l_n*56)+50], c3_9);
    __m128d c3_11 = _mm_load_sd(&C[(l_n*56)+52]);
    __m128d a3_11 = _mm_load_sd(&A[37]);
#if defined(__SSE3__) && defined(__AVX__)
    c3_11 = _mm_add_sd(c3_11, _mm_mul_sd(a3_11, _mm256_castpd256_pd128(b3)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c3_11 = _mm_add_sd(c3_11, _mm_mul_sd(a3_11, b3));
#endif
    _mm_store_sd(&C[(l_n*56)+52], c3_11);
#else
    C[(l_n*56)+4] += A[26] * B[(l_n*21)+3];
    C[(l_n*56)+5] += A[27] * B[(l_n*21)+3];
    C[(l_n*56)+6] += A[28] * B[(l_n*21)+3];
    C[(l_n*56)+14] += A[29] * B[(l_n*21)+3];
    C[(l_n*56)+15] += A[30] * B[(l_n*21)+3];
    C[(l_n*56)+16] += A[31] * B[(l_n*21)+3];
    C[(l_n*56)+29] += A[32] * B[(l_n*21)+3];
    C[(l_n*56)+30] += A[33] * B[(l_n*21)+3];
    C[(l_n*56)+31] += A[34] * B[(l_n*21)+3];
    C[(l_n*56)+50] += A[35] * B[(l_n*21)+3];
    C[(l_n*56)+51] += A[36] * B[(l_n*21)+3];
    C[(l_n*56)+52] += A[37] * B[(l_n*21)+3];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b4 = _mm256_broadcast_sd(&B[(l_n*21)+4]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b4 = _mm_loaddup_pd(&B[(l_n*21)+4]);
#endif
    __m128d c4_0 = _mm_loadu_pd(&C[(l_n*56)+4]);
    __m128d a4_0 = _mm_loadu_pd(&A[38]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_0 = _mm_add_pd(c4_0, _mm_mul_pd(a4_0, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_0 = _mm_add_pd(c4_0, _mm_mul_pd(a4_0, b4));
#endif
    _mm_storeu_pd(&C[(l_n*56)+4], c4_0);
    __m128d c4_2 = _mm_load_sd(&C[(l_n*56)+6]);
    __m128d a4_2 = _mm_load_sd(&A[40]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_2 = _mm_add_sd(c4_2, _mm_mul_sd(a4_2, b4));
#endif
    _mm_store_sd(&C[(l_n*56)+6], c4_2);
    __m128d c4_3 = _mm_loadu_pd(&C[(l_n*56)+14]);
    __m128d a4_3 = _mm_loadu_pd(&A[41]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_3 = _mm_add_pd(c4_3, _mm_mul_pd(a4_3, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_3 = _mm_add_pd(c4_3, _mm_mul_pd(a4_3, b4));
#endif
    _mm_storeu_pd(&C[(l_n*56)+14], c4_3);
    __m128d c4_5 = _mm_load_sd(&C[(l_n*56)+16]);
    __m128d a4_5 = _mm_load_sd(&A[43]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_5 = _mm_add_sd(c4_5, _mm_mul_sd(a4_5, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_5 = _mm_add_sd(c4_5, _mm_mul_sd(a4_5, b4));
#endif
    _mm_store_sd(&C[(l_n*56)+16], c4_5);
    __m128d c4_6 = _mm_loadu_pd(&C[(l_n*56)+29]);
    __m128d a4_6 = _mm_loadu_pd(&A[44]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_6 = _mm_add_pd(c4_6, _mm_mul_pd(a4_6, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_6 = _mm_add_pd(c4_6, _mm_mul_pd(a4_6, b4));
#endif
    _mm_storeu_pd(&C[(l_n*56)+29], c4_6);
    __m128d c4_8 = _mm_load_sd(&C[(l_n*56)+31]);
    __m128d a4_8 = _mm_load_sd(&A[46]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_8 = _mm_add_sd(c4_8, _mm_mul_sd(a4_8, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_8 = _mm_add_sd(c4_8, _mm_mul_sd(a4_8, b4));
#endif
    _mm_store_sd(&C[(l_n*56)+31], c4_8);
    __m128d c4_9 = _mm_loadu_pd(&C[(l_n*56)+50]);
    __m128d a4_9 = _mm_loadu_pd(&A[47]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_9 = _mm_add_pd(c4_9, _mm_mul_pd(a4_9, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_9 = _mm_add_pd(c4_9, _mm_mul_pd(a4_9, b4));
#endif
    _mm_storeu_pd(&C[(l_n*56)+50], c4_9);
    __m128d c4_11 = _mm_load_sd(&C[(l_n*56)+52]);
    __m128d a4_11 = _mm_load_sd(&A[49]);
#if defined(__SSE3__) && defined(__AVX__)
    c4_11 = _mm_add_sd(c4_11, _mm_mul_sd(a4_11, _mm256_castpd256_pd128(b4)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c4_11 = _mm_add_sd(c4_11, _mm_mul_sd(a4_11, b4));
#endif
    _mm_store_sd(&C[(l_n*56)+52], c4_11);
#else
    C[(l_n*56)+4] += A[38] * B[(l_n*21)+4];
    C[(l_n*56)+5] += A[39] * B[(l_n*21)+4];
    C[(l_n*56)+6] += A[40] * B[(l_n*21)+4];
    C[(l_n*56)+14] += A[41] * B[(l_n*21)+4];
    C[(l_n*56)+15] += A[42] * B[(l_n*21)+4];
    C[(l_n*56)+16] += A[43] * B[(l_n*21)+4];
    C[(l_n*56)+29] += A[44] * B[(l_n*21)+4];
    C[(l_n*56)+30] += A[45] * B[(l_n*21)+4];
    C[(l_n*56)+31] += A[46] * B[(l_n*21)+4];
    C[(l_n*56)+50] += A[47] * B[(l_n*21)+4];
    C[(l_n*56)+51] += A[48] * B[(l_n*21)+4];
    C[(l_n*56)+52] += A[49] * B[(l_n*21)+4];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b5 = _mm256_broadcast_sd(&B[(l_n*21)+5]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b5 = _mm_loaddup_pd(&B[(l_n*21)+5]);
#endif
    __m128d c5_0 = _mm_loadu_pd(&C[(l_n*56)+4]);
    __m128d a5_0 = _mm_loadu_pd(&A[50]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_0 = _mm_add_pd(c5_0, _mm_mul_pd(a5_0, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_0 = _mm_add_pd(c5_0, _mm_mul_pd(a5_0, b5));
#endif
    _mm_storeu_pd(&C[(l_n*56)+4], c5_0);
    __m128d c5_2 = _mm_load_sd(&C[(l_n*56)+6]);
    __m128d a5_2 = _mm_load_sd(&A[52]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_2 = _mm_add_sd(c5_2, _mm_mul_sd(a5_2, b5));
#endif
    _mm_store_sd(&C[(l_n*56)+6], c5_2);
    __m128d c5_3 = _mm_loadu_pd(&C[(l_n*56)+14]);
    __m128d a5_3 = _mm_loadu_pd(&A[53]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_3 = _mm_add_pd(c5_3, _mm_mul_pd(a5_3, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_3 = _mm_add_pd(c5_3, _mm_mul_pd(a5_3, b5));
#endif
    _mm_storeu_pd(&C[(l_n*56)+14], c5_3);
    __m128d c5_5 = _mm_load_sd(&C[(l_n*56)+16]);
    __m128d a5_5 = _mm_load_sd(&A[55]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_5 = _mm_add_sd(c5_5, _mm_mul_sd(a5_5, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_5 = _mm_add_sd(c5_5, _mm_mul_sd(a5_5, b5));
#endif
    _mm_store_sd(&C[(l_n*56)+16], c5_5);
    __m128d c5_6 = _mm_loadu_pd(&C[(l_n*56)+29]);
    __m128d a5_6 = _mm_loadu_pd(&A[56]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_6 = _mm_add_pd(c5_6, _mm_mul_pd(a5_6, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_6 = _mm_add_pd(c5_6, _mm_mul_pd(a5_6, b5));
#endif
    _mm_storeu_pd(&C[(l_n*56)+29], c5_6);
    __m128d c5_8 = _mm_load_sd(&C[(l_n*56)+31]);
    __m128d a5_8 = _mm_load_sd(&A[58]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_8 = _mm_add_sd(c5_8, _mm_mul_sd(a5_8, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_8 = _mm_add_sd(c5_8, _mm_mul_sd(a5_8, b5));
#endif
    _mm_store_sd(&C[(l_n*56)+31], c5_8);
    __m128d c5_9 = _mm_loadu_pd(&C[(l_n*56)+50]);
    __m128d a5_9 = _mm_loadu_pd(&A[59]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_9 = _mm_add_pd(c5_9, _mm_mul_pd(a5_9, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_9 = _mm_add_pd(c5_9, _mm_mul_pd(a5_9, b5));
#endif
    _mm_storeu_pd(&C[(l_n*56)+50], c5_9);
    __m128d c5_11 = _mm_load_sd(&C[(l_n*56)+52]);
    __m128d a5_11 = _mm_load_sd(&A[61]);
#if defined(__SSE3__) && defined(__AVX__)
    c5_11 = _mm_add_sd(c5_11, _mm_mul_sd(a5_11, _mm256_castpd256_pd128(b5)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c5_11 = _mm_add_sd(c5_11, _mm_mul_sd(a5_11, b5));
#endif
    _mm_store_sd(&C[(l_n*56)+52], c5_11);
#else
    C[(l_n*56)+4] += A[50] * B[(l_n*21)+5];
    C[(l_n*56)+5] += A[51] * B[(l_n*21)+5];
    C[(l_n*56)+6] += A[52] * B[(l_n*21)+5];
    C[(l_n*56)+14] += A[53] * B[(l_n*21)+5];
    C[(l_n*56)+15] += A[54] * B[(l_n*21)+5];
    C[(l_n*56)+16] += A[55] * B[(l_n*21)+5];
    C[(l_n*56)+29] += A[56] * B[(l_n*21)+5];
    C[(l_n*56)+30] += A[57] * B[(l_n*21)+5];
    C[(l_n*56)+31] += A[58] * B[(l_n*21)+5];
    C[(l_n*56)+50] += A[59] * B[(l_n*21)+5];
    C[(l_n*56)+51] += A[60] * B[(l_n*21)+5];
    C[(l_n*56)+52] += A[61] * B[(l_n*21)+5];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b6 = _mm256_broadcast_sd(&B[(l_n*21)+6]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b6 = _mm_loaddup_pd(&B[(l_n*21)+6]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c6_0 = _mm256_loadu_pd(&C[(l_n*56)+10]);
    __m256d a6_0 = _mm256_loadu_pd(&A[62]);
    c6_0 = _mm256_add_pd(c6_0, _mm256_mul_pd(a6_0, b6));
    _mm256_storeu_pd(&C[(l_n*56)+10], c6_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c6_0 = _mm_loadu_pd(&C[(l_n*56)+10]);
    __m128d a6_0 = _mm_loadu_pd(&A[62]);
    c6_0 = _mm_add_pd(c6_0, _mm_mul_pd(a6_0, b6));
    _mm_storeu_pd(&C[(l_n*56)+10], c6_0);
    __m128d c6_2 = _mm_loadu_pd(&C[(l_n*56)+12]);
    __m128d a6_2 = _mm_loadu_pd(&A[64]);
    c6_2 = _mm_add_pd(c6_2, _mm_mul_pd(a6_2, b6));
    _mm_storeu_pd(&C[(l_n*56)+12], c6_2);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c6_4 = _mm256_loadu_pd(&C[(l_n*56)+25]);
    __m256d a6_4 = _mm256_loadu_pd(&A[66]);
    c6_4 = _mm256_add_pd(c6_4, _mm256_mul_pd(a6_4, b6));
    _mm256_storeu_pd(&C[(l_n*56)+25], c6_4);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c6_4 = _mm_loadu_pd(&C[(l_n*56)+25]);
    __m128d a6_4 = _mm_loadu_pd(&A[66]);
    c6_4 = _mm_add_pd(c6_4, _mm_mul_pd(a6_4, b6));
    _mm_storeu_pd(&C[(l_n*56)+25], c6_4);
    __m128d c6_6 = _mm_loadu_pd(&C[(l_n*56)+27]);
    __m128d a6_6 = _mm_loadu_pd(&A[68]);
    c6_6 = _mm_add_pd(c6_6, _mm_mul_pd(a6_6, b6));
    _mm_storeu_pd(&C[(l_n*56)+27], c6_6);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c6_8 = _mm256_loadu_pd(&C[(l_n*56)+46]);
    __m256d a6_8 = _mm256_loadu_pd(&A[70]);
    c6_8 = _mm256_add_pd(c6_8, _mm256_mul_pd(a6_8, b6));
    _mm256_storeu_pd(&C[(l_n*56)+46], c6_8);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c6_8 = _mm_loadu_pd(&C[(l_n*56)+46]);
    __m128d a6_8 = _mm_loadu_pd(&A[70]);
    c6_8 = _mm_add_pd(c6_8, _mm_mul_pd(a6_8, b6));
    _mm_storeu_pd(&C[(l_n*56)+46], c6_8);
    __m128d c6_10 = _mm_loadu_pd(&C[(l_n*56)+48]);
    __m128d a6_10 = _mm_loadu_pd(&A[72]);
    c6_10 = _mm_add_pd(c6_10, _mm_mul_pd(a6_10, b6));
    _mm_storeu_pd(&C[(l_n*56)+48], c6_10);
#endif
#else
    C[(l_n*56)+10] += A[62] * B[(l_n*21)+6];
    C[(l_n*56)+11] += A[63] * B[(l_n*21)+6];
    C[(l_n*56)+12] += A[64] * B[(l_n*21)+6];
    C[(l_n*56)+13] += A[65] * B[(l_n*21)+6];
    C[(l_n*56)+25] += A[66] * B[(l_n*21)+6];
    C[(l_n*56)+26] += A[67] * B[(l_n*21)+6];
    C[(l_n*56)+27] += A[68] * B[(l_n*21)+6];
    C[(l_n*56)+28] += A[69] * B[(l_n*21)+6];
    C[(l_n*56)+46] += A[70] * B[(l_n*21)+6];
    C[(l_n*56)+47] += A[71] * B[(l_n*21)+6];
    C[(l_n*56)+48] += A[72] * B[(l_n*21)+6];
    C[(l_n*56)+49] += A[73] * B[(l_n*21)+6];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b7 = _mm256_broadcast_sd(&B[(l_n*21)+7]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b7 = _mm_loaddup_pd(&B[(l_n*21)+7]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c7_0 = _mm256_loadu_pd(&C[(l_n*56)+10]);
    __m256d a7_0 = _mm256_loadu_pd(&A[74]);
    c7_0 = _mm256_add_pd(c7_0, _mm256_mul_pd(a7_0, b7));
    _mm256_storeu_pd(&C[(l_n*56)+10], c7_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c7_0 = _mm_loadu_pd(&C[(l_n*56)+10]);
    __m128d a7_0 = _mm_loadu_pd(&A[74]);
    c7_0 = _mm_add_pd(c7_0, _mm_mul_pd(a7_0, b7));
    _mm_storeu_pd(&C[(l_n*56)+10], c7_0);
    __m128d c7_2 = _mm_loadu_pd(&C[(l_n*56)+12]);
    __m128d a7_2 = _mm_loadu_pd(&A[76]);
    c7_2 = _mm_add_pd(c7_2, _mm_mul_pd(a7_2, b7));
    _mm_storeu_pd(&C[(l_n*56)+12], c7_2);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c7_4 = _mm256_loadu_pd(&C[(l_n*56)+25]);
    __m256d a7_4 = _mm256_loadu_pd(&A[78]);
    c7_4 = _mm256_add_pd(c7_4, _mm256_mul_pd(a7_4, b7));
    _mm256_storeu_pd(&C[(l_n*56)+25], c7_4);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c7_4 = _mm_loadu_pd(&C[(l_n*56)+25]);
    __m128d a7_4 = _mm_loadu_pd(&A[78]);
    c7_4 = _mm_add_pd(c7_4, _mm_mul_pd(a7_4, b7));
    _mm_storeu_pd(&C[(l_n*56)+25], c7_4);
    __m128d c7_6 = _mm_loadu_pd(&C[(l_n*56)+27]);
    __m128d a7_6 = _mm_loadu_pd(&A[80]);
    c7_6 = _mm_add_pd(c7_6, _mm_mul_pd(a7_6, b7));
    _mm_storeu_pd(&C[(l_n*56)+27], c7_6);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c7_8 = _mm256_loadu_pd(&C[(l_n*56)+46]);
    __m256d a7_8 = _mm256_loadu_pd(&A[82]);
    c7_8 = _mm256_add_pd(c7_8, _mm256_mul_pd(a7_8, b7));
    _mm256_storeu_pd(&C[(l_n*56)+46], c7_8);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c7_8 = _mm_loadu_pd(&C[(l_n*56)+46]);
    __m128d a7_8 = _mm_loadu_pd(&A[82]);
    c7_8 = _mm_add_pd(c7_8, _mm_mul_pd(a7_8, b7));
    _mm_storeu_pd(&C[(l_n*56)+46], c7_8);
    __m128d c7_10 = _mm_loadu_pd(&C[(l_n*56)+48]);
    __m128d a7_10 = _mm_loadu_pd(&A[84]);
    c7_10 = _mm_add_pd(c7_10, _mm_mul_pd(a7_10, b7));
    _mm_storeu_pd(&C[(l_n*56)+48], c7_10);
#endif
#else
    C[(l_n*56)+10] += A[74] * B[(l_n*21)+7];
    C[(l_n*56)+11] += A[75] * B[(l_n*21)+7];
    C[(l_n*56)+12] += A[76] * B[(l_n*21)+7];
    C[(l_n*56)+13] += A[77] * B[(l_n*21)+7];
    C[(l_n*56)+25] += A[78] * B[(l_n*21)+7];
    C[(l_n*56)+26] += A[79] * B[(l_n*21)+7];
    C[(l_n*56)+27] += A[80] * B[(l_n*21)+7];
    C[(l_n*56)+28] += A[81] * B[(l_n*21)+7];
    C[(l_n*56)+46] += A[82] * B[(l_n*21)+7];
    C[(l_n*56)+47] += A[83] * B[(l_n*21)+7];
    C[(l_n*56)+48] += A[84] * B[(l_n*21)+7];
    C[(l_n*56)+49] += A[85] * B[(l_n*21)+7];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b8 = _mm256_broadcast_sd(&B[(l_n*21)+8]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b8 = _mm_loaddup_pd(&B[(l_n*21)+8]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c8_0 = _mm256_loadu_pd(&C[(l_n*56)+10]);
    __m256d a8_0 = _mm256_loadu_pd(&A[86]);
    c8_0 = _mm256_add_pd(c8_0, _mm256_mul_pd(a8_0, b8));
    _mm256_storeu_pd(&C[(l_n*56)+10], c8_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c8_0 = _mm_loadu_pd(&C[(l_n*56)+10]);
    __m128d a8_0 = _mm_loadu_pd(&A[86]);
    c8_0 = _mm_add_pd(c8_0, _mm_mul_pd(a8_0, b8));
    _mm_storeu_pd(&C[(l_n*56)+10], c8_0);
    __m128d c8_2 = _mm_loadu_pd(&C[(l_n*56)+12]);
    __m128d a8_2 = _mm_loadu_pd(&A[88]);
    c8_2 = _mm_add_pd(c8_2, _mm_mul_pd(a8_2, b8));
    _mm_storeu_pd(&C[(l_n*56)+12], c8_2);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c8_4 = _mm256_loadu_pd(&C[(l_n*56)+25]);
    __m256d a8_4 = _mm256_loadu_pd(&A[90]);
    c8_4 = _mm256_add_pd(c8_4, _mm256_mul_pd(a8_4, b8));
    _mm256_storeu_pd(&C[(l_n*56)+25], c8_4);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c8_4 = _mm_loadu_pd(&C[(l_n*56)+25]);
    __m128d a8_4 = _mm_loadu_pd(&A[90]);
    c8_4 = _mm_add_pd(c8_4, _mm_mul_pd(a8_4, b8));
    _mm_storeu_pd(&C[(l_n*56)+25], c8_4);
    __m128d c8_6 = _mm_loadu_pd(&C[(l_n*56)+27]);
    __m128d a8_6 = _mm_loadu_pd(&A[92]);
    c8_6 = _mm_add_pd(c8_6, _mm_mul_pd(a8_6, b8));
    _mm_storeu_pd(&C[(l_n*56)+27], c8_6);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c8_8 = _mm256_loadu_pd(&C[(l_n*56)+46]);
    __m256d a8_8 = _mm256_loadu_pd(&A[94]);
    c8_8 = _mm256_add_pd(c8_8, _mm256_mul_pd(a8_8, b8));
    _mm256_storeu_pd(&C[(l_n*56)+46], c8_8);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c8_8 = _mm_loadu_pd(&C[(l_n*56)+46]);
    __m128d a8_8 = _mm_loadu_pd(&A[94]);
    c8_8 = _mm_add_pd(c8_8, _mm_mul_pd(a8_8, b8));
    _mm_storeu_pd(&C[(l_n*56)+46], c8_8);
    __m128d c8_10 = _mm_loadu_pd(&C[(l_n*56)+48]);
    __m128d a8_10 = _mm_loadu_pd(&A[96]);
    c8_10 = _mm_add_pd(c8_10, _mm_mul_pd(a8_10, b8));
    _mm_storeu_pd(&C[(l_n*56)+48], c8_10);
#endif
#else
    C[(l_n*56)+10] += A[86] * B[(l_n*21)+8];
    C[(l_n*56)+11] += A[87] * B[(l_n*21)+8];
    C[(l_n*56)+12] += A[88] * B[(l_n*21)+8];
    C[(l_n*56)+13] += A[89] * B[(l_n*21)+8];
    C[(l_n*56)+25] += A[90] * B[(l_n*21)+8];
    C[(l_n*56)+26] += A[91] * B[(l_n*21)+8];
    C[(l_n*56)+27] += A[92] * B[(l_n*21)+8];
    C[(l_n*56)+28] += A[93] * B[(l_n*21)+8];
    C[(l_n*56)+46] += A[94] * B[(l_n*21)+8];
    C[(l_n*56)+47] += A[95] * B[(l_n*21)+8];
    C[(l_n*56)+48] += A[96] * B[(l_n*21)+8];
    C[(l_n*56)+49] += A[97] * B[(l_n*21)+8];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b9 = _mm256_broadcast_sd(&B[(l_n*21)+9]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b9 = _mm_loaddup_pd(&B[(l_n*21)+9]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c9_0 = _mm256_loadu_pd(&C[(l_n*56)+10]);
    __m256d a9_0 = _mm256_loadu_pd(&A[98]);
    c9_0 = _mm256_add_pd(c9_0, _mm256_mul_pd(a9_0, b9));
    _mm256_storeu_pd(&C[(l_n*56)+10], c9_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c9_0 = _mm_loadu_pd(&C[(l_n*56)+10]);
    __m128d a9_0 = _mm_loadu_pd(&A[98]);
    c9_0 = _mm_add_pd(c9_0, _mm_mul_pd(a9_0, b9));
    _mm_storeu_pd(&C[(l_n*56)+10], c9_0);
    __m128d c9_2 = _mm_loadu_pd(&C[(l_n*56)+12]);
    __m128d a9_2 = _mm_loadu_pd(&A[100]);
    c9_2 = _mm_add_pd(c9_2, _mm_mul_pd(a9_2, b9));
    _mm_storeu_pd(&C[(l_n*56)+12], c9_2);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c9_4 = _mm256_loadu_pd(&C[(l_n*56)+25]);
    __m256d a9_4 = _mm256_loadu_pd(&A[102]);
    c9_4 = _mm256_add_pd(c9_4, _mm256_mul_pd(a9_4, b9));
    _mm256_storeu_pd(&C[(l_n*56)+25], c9_4);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c9_4 = _mm_loadu_pd(&C[(l_n*56)+25]);
    __m128d a9_4 = _mm_loadu_pd(&A[102]);
    c9_4 = _mm_add_pd(c9_4, _mm_mul_pd(a9_4, b9));
    _mm_storeu_pd(&C[(l_n*56)+25], c9_4);
    __m128d c9_6 = _mm_loadu_pd(&C[(l_n*56)+27]);
    __m128d a9_6 = _mm_loadu_pd(&A[104]);
    c9_6 = _mm_add_pd(c9_6, _mm_mul_pd(a9_6, b9));
    _mm_storeu_pd(&C[(l_n*56)+27], c9_6);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c9_8 = _mm256_loadu_pd(&C[(l_n*56)+46]);
    __m256d a9_8 = _mm256_loadu_pd(&A[106]);
    c9_8 = _mm256_add_pd(c9_8, _mm256_mul_pd(a9_8, b9));
    _mm256_storeu_pd(&C[(l_n*56)+46], c9_8);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c9_8 = _mm_loadu_pd(&C[(l_n*56)+46]);
    __m128d a9_8 = _mm_loadu_pd(&A[106]);
    c9_8 = _mm_add_pd(c9_8, _mm_mul_pd(a9_8, b9));
    _mm_storeu_pd(&C[(l_n*56)+46], c9_8);
    __m128d c9_10 = _mm_loadu_pd(&C[(l_n*56)+48]);
    __m128d a9_10 = _mm_loadu_pd(&A[108]);
    c9_10 = _mm_add_pd(c9_10, _mm_mul_pd(a9_10, b9));
    _mm_storeu_pd(&C[(l_n*56)+48], c9_10);
#endif
#else
    C[(l_n*56)+10] += A[98] * B[(l_n*21)+9];
    C[(l_n*56)+11] += A[99] * B[(l_n*21)+9];
    C[(l_n*56)+12] += A[100] * B[(l_n*21)+9];
    C[(l_n*56)+13] += A[101] * B[(l_n*21)+9];
    C[(l_n*56)+25] += A[102] * B[(l_n*21)+9];
    C[(l_n*56)+26] += A[103] * B[(l_n*21)+9];
    C[(l_n*56)+27] += A[104] * B[(l_n*21)+9];
    C[(l_n*56)+28] += A[105] * B[(l_n*21)+9];
    C[(l_n*56)+46] += A[106] * B[(l_n*21)+9];
    C[(l_n*56)+47] += A[107] * B[(l_n*21)+9];
    C[(l_n*56)+48] += A[108] * B[(l_n*21)+9];
    C[(l_n*56)+49] += A[109] * B[(l_n*21)+9];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b10 = _mm256_broadcast_sd(&B[(l_n*21)+10]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b10 = _mm_loaddup_pd(&B[(l_n*21)+10]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c10_0 = _mm256_loadu_pd(&C[(l_n*56)+20]);
    __m256d a10_0 = _mm256_loadu_pd(&A[110]);
    c10_0 = _mm256_add_pd(c10_0, _mm256_mul_pd(a10_0, b10));
    _mm256_storeu_pd(&C[(l_n*56)+20], c10_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c10_0 = _mm_loadu_pd(&C[(l_n*56)+20]);
    __m128d a10_0 = _mm_loadu_pd(&A[110]);
    c10_0 = _mm_add_pd(c10_0, _mm_mul_pd(a10_0, b10));
    _mm_storeu_pd(&C[(l_n*56)+20], c10_0);
    __m128d c10_2 = _mm_loadu_pd(&C[(l_n*56)+22]);
    __m128d a10_2 = _mm_loadu_pd(&A[112]);
    c10_2 = _mm_add_pd(c10_2, _mm_mul_pd(a10_2, b10));
    _mm_storeu_pd(&C[(l_n*56)+22], c10_2);
#endif
    __m128d c10_4 = _mm_load_sd(&C[(l_n*56)+24]);
    __m128d a10_4 = _mm_load_sd(&A[114]);
#if defined(__SSE3__) && defined(__AVX__)
    c10_4 = _mm_add_sd(c10_4, _mm_mul_sd(a10_4, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c10_4 = _mm_add_sd(c10_4, _mm_mul_sd(a10_4, b10));
#endif
    _mm_store_sd(&C[(l_n*56)+24], c10_4);
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c10_5 = _mm256_loadu_pd(&C[(l_n*56)+41]);
    __m256d a10_5 = _mm256_loadu_pd(&A[115]);
    c10_5 = _mm256_add_pd(c10_5, _mm256_mul_pd(a10_5, b10));
    _mm256_storeu_pd(&C[(l_n*56)+41], c10_5);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c10_5 = _mm_loadu_pd(&C[(l_n*56)+41]);
    __m128d a10_5 = _mm_loadu_pd(&A[115]);
    c10_5 = _mm_add_pd(c10_5, _mm_mul_pd(a10_5, b10));
    _mm_storeu_pd(&C[(l_n*56)+41], c10_5);
    __m128d c10_7 = _mm_loadu_pd(&C[(l_n*56)+43]);
    __m128d a10_7 = _mm_loadu_pd(&A[117]);
    c10_7 = _mm_add_pd(c10_7, _mm_mul_pd(a10_7, b10));
    _mm_storeu_pd(&C[(l_n*56)+43], c10_7);
#endif
    __m128d c10_9 = _mm_load_sd(&C[(l_n*56)+45]);
    __m128d a10_9 = _mm_load_sd(&A[119]);
#if defined(__SSE3__) && defined(__AVX__)
    c10_9 = _mm_add_sd(c10_9, _mm_mul_sd(a10_9, _mm256_castpd256_pd128(b10)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c10_9 = _mm_add_sd(c10_9, _mm_mul_sd(a10_9, b10));
#endif
    _mm_store_sd(&C[(l_n*56)+45], c10_9);
#else
    C[(l_n*56)+20] += A[110] * B[(l_n*21)+10];
    C[(l_n*56)+21] += A[111] * B[(l_n*21)+10];
    C[(l_n*56)+22] += A[112] * B[(l_n*21)+10];
    C[(l_n*56)+23] += A[113] * B[(l_n*21)+10];
    C[(l_n*56)+24] += A[114] * B[(l_n*21)+10];
    C[(l_n*56)+41] += A[115] * B[(l_n*21)+10];
    C[(l_n*56)+42] += A[116] * B[(l_n*21)+10];
    C[(l_n*56)+43] += A[117] * B[(l_n*21)+10];
    C[(l_n*56)+44] += A[118] * B[(l_n*21)+10];
    C[(l_n*56)+45] += A[119] * B[(l_n*21)+10];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b11 = _mm256_broadcast_sd(&B[(l_n*21)+11]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b11 = _mm_loaddup_pd(&B[(l_n*21)+11]);
#endif
    __m128d c11_0 = _mm_loadu_pd(&C[(l_n*56)+20]);
    __m128d a11_0 = _mm_loadu_pd(&A[120]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_0 = _mm_add_pd(c11_0, _mm_mul_pd(a11_0, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_0 = _mm_add_pd(c11_0, _mm_mul_pd(a11_0, b11));
#endif
    _mm_storeu_pd(&C[(l_n*56)+20], c11_0);
    __m128d c11_2 = _mm_load_sd(&C[(l_n*56)+22]);
    __m128d a11_2 = _mm_load_sd(&A[122]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_2 = _mm_add_sd(c11_2, _mm_mul_sd(a11_2, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_2 = _mm_add_sd(c11_2, _mm_mul_sd(a11_2, b11));
#endif
    _mm_store_sd(&C[(l_n*56)+22], c11_2);
    __m128d c11_3 = _mm_load_sd(&C[(l_n*56)+24]);
    __m128d a11_3 = _mm_load_sd(&A[123]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_3 = _mm_add_sd(c11_3, _mm_mul_sd(a11_3, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_3 = _mm_add_sd(c11_3, _mm_mul_sd(a11_3, b11));
#endif
    _mm_store_sd(&C[(l_n*56)+24], c11_3);
    __m128d c11_4 = _mm_loadu_pd(&C[(l_n*56)+41]);
    __m128d a11_4 = _mm_loadu_pd(&A[124]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_4 = _mm_add_pd(c11_4, _mm_mul_pd(a11_4, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_4 = _mm_add_pd(c11_4, _mm_mul_pd(a11_4, b11));
#endif
    _mm_storeu_pd(&C[(l_n*56)+41], c11_4);
    __m128d c11_6 = _mm_load_sd(&C[(l_n*56)+43]);
    __m128d a11_6 = _mm_load_sd(&A[126]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_6 = _mm_add_sd(c11_6, _mm_mul_sd(a11_6, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_6 = _mm_add_sd(c11_6, _mm_mul_sd(a11_6, b11));
#endif
    _mm_store_sd(&C[(l_n*56)+43], c11_6);
    __m128d c11_7 = _mm_load_sd(&C[(l_n*56)+45]);
    __m128d a11_7 = _mm_load_sd(&A[127]);
#if defined(__SSE3__) && defined(__AVX__)
    c11_7 = _mm_add_sd(c11_7, _mm_mul_sd(a11_7, _mm256_castpd256_pd128(b11)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c11_7 = _mm_add_sd(c11_7, _mm_mul_sd(a11_7, b11));
#endif
    _mm_store_sd(&C[(l_n*56)+45], c11_7);
#else
    C[(l_n*56)+20] += A[120] * B[(l_n*21)+11];
    C[(l_n*56)+21] += A[121] * B[(l_n*21)+11];
    C[(l_n*56)+22] += A[122] * B[(l_n*21)+11];
    C[(l_n*56)+24] += A[123] * B[(l_n*21)+11];
    C[(l_n*56)+41] += A[124] * B[(l_n*21)+11];
    C[(l_n*56)+42] += A[125] * B[(l_n*21)+11];
    C[(l_n*56)+43] += A[126] * B[(l_n*21)+11];
    C[(l_n*56)+45] += A[127] * B[(l_n*21)+11];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b12 = _mm256_broadcast_sd(&B[(l_n*21)+12]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b12 = _mm_loaddup_pd(&B[(l_n*21)+12]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c12_0 = _mm256_loadu_pd(&C[(l_n*56)+20]);
    __m256d a12_0 = _mm256_loadu_pd(&A[128]);
    c12_0 = _mm256_add_pd(c12_0, _mm256_mul_pd(a12_0, b12));
    _mm256_storeu_pd(&C[(l_n*56)+20], c12_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c12_0 = _mm_loadu_pd(&C[(l_n*56)+20]);
    __m128d a12_0 = _mm_loadu_pd(&A[128]);
    c12_0 = _mm_add_pd(c12_0, _mm_mul_pd(a12_0, b12));
    _mm_storeu_pd(&C[(l_n*56)+20], c12_0);
    __m128d c12_2 = _mm_loadu_pd(&C[(l_n*56)+22]);
    __m128d a12_2 = _mm_loadu_pd(&A[130]);
    c12_2 = _mm_add_pd(c12_2, _mm_mul_pd(a12_2, b12));
    _mm_storeu_pd(&C[(l_n*56)+22], c12_2);
#endif
    __m128d c12_4 = _mm_load_sd(&C[(l_n*56)+24]);
    __m128d a12_4 = _mm_load_sd(&A[132]);
#if defined(__SSE3__) && defined(__AVX__)
    c12_4 = _mm_add_sd(c12_4, _mm_mul_sd(a12_4, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c12_4 = _mm_add_sd(c12_4, _mm_mul_sd(a12_4, b12));
#endif
    _mm_store_sd(&C[(l_n*56)+24], c12_4);
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c12_5 = _mm256_loadu_pd(&C[(l_n*56)+41]);
    __m256d a12_5 = _mm256_loadu_pd(&A[133]);
    c12_5 = _mm256_add_pd(c12_5, _mm256_mul_pd(a12_5, b12));
    _mm256_storeu_pd(&C[(l_n*56)+41], c12_5);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c12_5 = _mm_loadu_pd(&C[(l_n*56)+41]);
    __m128d a12_5 = _mm_loadu_pd(&A[133]);
    c12_5 = _mm_add_pd(c12_5, _mm_mul_pd(a12_5, b12));
    _mm_storeu_pd(&C[(l_n*56)+41], c12_5);
    __m128d c12_7 = _mm_loadu_pd(&C[(l_n*56)+43]);
    __m128d a12_7 = _mm_loadu_pd(&A[135]);
    c12_7 = _mm_add_pd(c12_7, _mm_mul_pd(a12_7, b12));
    _mm_storeu_pd(&C[(l_n*56)+43], c12_7);
#endif
    __m128d c12_9 = _mm_load_sd(&C[(l_n*56)+45]);
    __m128d a12_9 = _mm_load_sd(&A[137]);
#if defined(__SSE3__) && defined(__AVX__)
    c12_9 = _mm_add_sd(c12_9, _mm_mul_sd(a12_9, _mm256_castpd256_pd128(b12)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c12_9 = _mm_add_sd(c12_9, _mm_mul_sd(a12_9, b12));
#endif
    _mm_store_sd(&C[(l_n*56)+45], c12_9);
#else
    C[(l_n*56)+20] += A[128] * B[(l_n*21)+12];
    C[(l_n*56)+21] += A[129] * B[(l_n*21)+12];
    C[(l_n*56)+22] += A[130] * B[(l_n*21)+12];
    C[(l_n*56)+23] += A[131] * B[(l_n*21)+12];
    C[(l_n*56)+24] += A[132] * B[(l_n*21)+12];
    C[(l_n*56)+41] += A[133] * B[(l_n*21)+12];
    C[(l_n*56)+42] += A[134] * B[(l_n*21)+12];
    C[(l_n*56)+43] += A[135] * B[(l_n*21)+12];
    C[(l_n*56)+44] += A[136] * B[(l_n*21)+12];
    C[(l_n*56)+45] += A[137] * B[(l_n*21)+12];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b13 = _mm256_broadcast_sd(&B[(l_n*21)+13]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b13 = _mm_loaddup_pd(&B[(l_n*21)+13]);
#endif
    __m128d c13_0 = _mm_load_sd(&C[(l_n*56)+20]);
    __m128d a13_0 = _mm_load_sd(&A[138]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_0 = _mm_add_sd(c13_0, _mm_mul_sd(a13_0, b13));
#endif
    _mm_store_sd(&C[(l_n*56)+20], c13_0);
    __m128d c13_1 = _mm_loadu_pd(&C[(l_n*56)+22]);
    __m128d a13_1 = _mm_loadu_pd(&A[139]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_1 = _mm_add_pd(c13_1, _mm_mul_pd(a13_1, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_1 = _mm_add_pd(c13_1, _mm_mul_pd(a13_1, b13));
#endif
    _mm_storeu_pd(&C[(l_n*56)+22], c13_1);
    __m128d c13_3 = _mm_load_sd(&C[(l_n*56)+24]);
    __m128d a13_3 = _mm_load_sd(&A[141]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_3 = _mm_add_sd(c13_3, _mm_mul_sd(a13_3, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_3 = _mm_add_sd(c13_3, _mm_mul_sd(a13_3, b13));
#endif
    _mm_store_sd(&C[(l_n*56)+24], c13_3);
    __m128d c13_4 = _mm_load_sd(&C[(l_n*56)+41]);
    __m128d a13_4 = _mm_load_sd(&A[142]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_4 = _mm_add_sd(c13_4, _mm_mul_sd(a13_4, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_4 = _mm_add_sd(c13_4, _mm_mul_sd(a13_4, b13));
#endif
    _mm_store_sd(&C[(l_n*56)+41], c13_4);
    __m128d c13_5 = _mm_loadu_pd(&C[(l_n*56)+43]);
    __m128d a13_5 = _mm_loadu_pd(&A[143]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_5 = _mm_add_pd(c13_5, _mm_mul_pd(a13_5, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_5 = _mm_add_pd(c13_5, _mm_mul_pd(a13_5, b13));
#endif
    _mm_storeu_pd(&C[(l_n*56)+43], c13_5);
    __m128d c13_7 = _mm_load_sd(&C[(l_n*56)+45]);
    __m128d a13_7 = _mm_load_sd(&A[145]);
#if defined(__SSE3__) && defined(__AVX__)
    c13_7 = _mm_add_sd(c13_7, _mm_mul_sd(a13_7, _mm256_castpd256_pd128(b13)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c13_7 = _mm_add_sd(c13_7, _mm_mul_sd(a13_7, b13));
#endif
    _mm_store_sd(&C[(l_n*56)+45], c13_7);
#else
    C[(l_n*56)+20] += A[138] * B[(l_n*21)+13];
    C[(l_n*56)+22] += A[139] * B[(l_n*21)+13];
    C[(l_n*56)+23] += A[140] * B[(l_n*21)+13];
    C[(l_n*56)+24] += A[141] * B[(l_n*21)+13];
    C[(l_n*56)+41] += A[142] * B[(l_n*21)+13];
    C[(l_n*56)+43] += A[143] * B[(l_n*21)+13];
    C[(l_n*56)+44] += A[144] * B[(l_n*21)+13];
    C[(l_n*56)+45] += A[145] * B[(l_n*21)+13];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b14 = _mm256_broadcast_sd(&B[(l_n*21)+14]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b14 = _mm_loaddup_pd(&B[(l_n*21)+14]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c14_0 = _mm256_loadu_pd(&C[(l_n*56)+20]);
    __m256d a14_0 = _mm256_loadu_pd(&A[146]);
    c14_0 = _mm256_add_pd(c14_0, _mm256_mul_pd(a14_0, b14));
    _mm256_storeu_pd(&C[(l_n*56)+20], c14_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c14_0 = _mm_loadu_pd(&C[(l_n*56)+20]);
    __m128d a14_0 = _mm_loadu_pd(&A[146]);
    c14_0 = _mm_add_pd(c14_0, _mm_mul_pd(a14_0, b14));
    _mm_storeu_pd(&C[(l_n*56)+20], c14_0);
    __m128d c14_2 = _mm_loadu_pd(&C[(l_n*56)+22]);
    __m128d a14_2 = _mm_loadu_pd(&A[148]);
    c14_2 = _mm_add_pd(c14_2, _mm_mul_pd(a14_2, b14));
    _mm_storeu_pd(&C[(l_n*56)+22], c14_2);
#endif
    __m128d c14_4 = _mm_load_sd(&C[(l_n*56)+24]);
    __m128d a14_4 = _mm_load_sd(&A[150]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_4 = _mm_add_sd(c14_4, _mm_mul_sd(a14_4, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_4 = _mm_add_sd(c14_4, _mm_mul_sd(a14_4, b14));
#endif
    _mm_store_sd(&C[(l_n*56)+24], c14_4);
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c14_5 = _mm256_loadu_pd(&C[(l_n*56)+41]);
    __m256d a14_5 = _mm256_loadu_pd(&A[151]);
    c14_5 = _mm256_add_pd(c14_5, _mm256_mul_pd(a14_5, b14));
    _mm256_storeu_pd(&C[(l_n*56)+41], c14_5);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c14_5 = _mm_loadu_pd(&C[(l_n*56)+41]);
    __m128d a14_5 = _mm_loadu_pd(&A[151]);
    c14_5 = _mm_add_pd(c14_5, _mm_mul_pd(a14_5, b14));
    _mm_storeu_pd(&C[(l_n*56)+41], c14_5);
    __m128d c14_7 = _mm_loadu_pd(&C[(l_n*56)+43]);
    __m128d a14_7 = _mm_loadu_pd(&A[153]);
    c14_7 = _mm_add_pd(c14_7, _mm_mul_pd(a14_7, b14));
    _mm_storeu_pd(&C[(l_n*56)+43], c14_7);
#endif
    __m128d c14_9 = _mm_load_sd(&C[(l_n*56)+45]);
    __m128d a14_9 = _mm_load_sd(&A[155]);
#if defined(__SSE3__) && defined(__AVX__)
    c14_9 = _mm_add_sd(c14_9, _mm_mul_sd(a14_9, _mm256_castpd256_pd128(b14)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c14_9 = _mm_add_sd(c14_9, _mm_mul_sd(a14_9, b14));
#endif
    _mm_store_sd(&C[(l_n*56)+45], c14_9);
#else
    C[(l_n*56)+20] += A[146] * B[(l_n*21)+14];
    C[(l_n*56)+21] += A[147] * B[(l_n*21)+14];
    C[(l_n*56)+22] += A[148] * B[(l_n*21)+14];
    C[(l_n*56)+23] += A[149] * B[(l_n*21)+14];
    C[(l_n*56)+24] += A[150] * B[(l_n*21)+14];
    C[(l_n*56)+41] += A[151] * B[(l_n*21)+14];
    C[(l_n*56)+42] += A[152] * B[(l_n*21)+14];
    C[(l_n*56)+43] += A[153] * B[(l_n*21)+14];
    C[(l_n*56)+44] += A[154] * B[(l_n*21)+14];
    C[(l_n*56)+45] += A[155] * B[(l_n*21)+14];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b15 = _mm256_broadcast_sd(&B[(l_n*21)+15]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b15 = _mm_loaddup_pd(&B[(l_n*21)+15]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c15_0 = _mm256_loadu_pd(&C[(l_n*56)+35]);
    __m256d a15_0 = _mm256_loadu_pd(&A[156]);
    c15_0 = _mm256_add_pd(c15_0, _mm256_mul_pd(a15_0, b15));
    _mm256_storeu_pd(&C[(l_n*56)+35], c15_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c15_0 = _mm_loadu_pd(&C[(l_n*56)+35]);
    __m128d a15_0 = _mm_loadu_pd(&A[156]);
    c15_0 = _mm_add_pd(c15_0, _mm_mul_pd(a15_0, b15));
    _mm_storeu_pd(&C[(l_n*56)+35], c15_0);
    __m128d c15_2 = _mm_loadu_pd(&C[(l_n*56)+37]);
    __m128d a15_2 = _mm_loadu_pd(&A[158]);
    c15_2 = _mm_add_pd(c15_2, _mm_mul_pd(a15_2, b15));
    _mm_storeu_pd(&C[(l_n*56)+37], c15_2);
#endif
    __m128d c15_4 = _mm_loadu_pd(&C[(l_n*56)+39]);
    __m128d a15_4 = _mm_loadu_pd(&A[160]);
#if defined(__SSE3__) && defined(__AVX__)
    c15_4 = _mm_add_pd(c15_4, _mm_mul_pd(a15_4, _mm256_castpd256_pd128(b15)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c15_4 = _mm_add_pd(c15_4, _mm_mul_pd(a15_4, b15));
#endif
    _mm_storeu_pd(&C[(l_n*56)+39], c15_4);
#else
    C[(l_n*56)+35] += A[156] * B[(l_n*21)+15];
    C[(l_n*56)+36] += A[157] * B[(l_n*21)+15];
    C[(l_n*56)+37] += A[158] * B[(l_n*21)+15];
    C[(l_n*56)+38] += A[159] * B[(l_n*21)+15];
    C[(l_n*56)+39] += A[160] * B[(l_n*21)+15];
    C[(l_n*56)+40] += A[161] * B[(l_n*21)+15];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b16 = _mm256_broadcast_sd(&B[(l_n*21)+16]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b16 = _mm_loaddup_pd(&B[(l_n*21)+16]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c16_0 = _mm256_loadu_pd(&C[(l_n*56)+35]);
    __m256d a16_0 = _mm256_loadu_pd(&A[162]);
    c16_0 = _mm256_add_pd(c16_0, _mm256_mul_pd(a16_0, b16));
    _mm256_storeu_pd(&C[(l_n*56)+35], c16_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c16_0 = _mm_loadu_pd(&C[(l_n*56)+35]);
    __m128d a16_0 = _mm_loadu_pd(&A[162]);
    c16_0 = _mm_add_pd(c16_0, _mm_mul_pd(a16_0, b16));
    _mm_storeu_pd(&C[(l_n*56)+35], c16_0);
    __m128d c16_2 = _mm_loadu_pd(&C[(l_n*56)+37]);
    __m128d a16_2 = _mm_loadu_pd(&A[164]);
    c16_2 = _mm_add_pd(c16_2, _mm_mul_pd(a16_2, b16));
    _mm_storeu_pd(&C[(l_n*56)+37], c16_2);
#endif
    __m128d c16_4 = _mm_loadu_pd(&C[(l_n*56)+39]);
    __m128d a16_4 = _mm_loadu_pd(&A[166]);
#if defined(__SSE3__) && defined(__AVX__)
    c16_4 = _mm_add_pd(c16_4, _mm_mul_pd(a16_4, _mm256_castpd256_pd128(b16)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c16_4 = _mm_add_pd(c16_4, _mm_mul_pd(a16_4, b16));
#endif
    _mm_storeu_pd(&C[(l_n*56)+39], c16_4);
#else
    C[(l_n*56)+35] += A[162] * B[(l_n*21)+16];
    C[(l_n*56)+36] += A[163] * B[(l_n*21)+16];
    C[(l_n*56)+37] += A[164] * B[(l_n*21)+16];
    C[(l_n*56)+38] += A[165] * B[(l_n*21)+16];
    C[(l_n*56)+39] += A[166] * B[(l_n*21)+16];
    C[(l_n*56)+40] += A[167] * B[(l_n*21)+16];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b17 = _mm256_broadcast_sd(&B[(l_n*21)+17]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b17 = _mm_loaddup_pd(&B[(l_n*21)+17]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c17_0 = _mm256_loadu_pd(&C[(l_n*56)+35]);
    __m256d a17_0 = _mm256_loadu_pd(&A[168]);
    c17_0 = _mm256_add_pd(c17_0, _mm256_mul_pd(a17_0, b17));
    _mm256_storeu_pd(&C[(l_n*56)+35], c17_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c17_0 = _mm_loadu_pd(&C[(l_n*56)+35]);
    __m128d a17_0 = _mm_loadu_pd(&A[168]);
    c17_0 = _mm_add_pd(c17_0, _mm_mul_pd(a17_0, b17));
    _mm_storeu_pd(&C[(l_n*56)+35], c17_0);
    __m128d c17_2 = _mm_loadu_pd(&C[(l_n*56)+37]);
    __m128d a17_2 = _mm_loadu_pd(&A[170]);
    c17_2 = _mm_add_pd(c17_2, _mm_mul_pd(a17_2, b17));
    _mm_storeu_pd(&C[(l_n*56)+37], c17_2);
#endif
    __m128d c17_4 = _mm_loadu_pd(&C[(l_n*56)+39]);
    __m128d a17_4 = _mm_loadu_pd(&A[172]);
#if defined(__SSE3__) && defined(__AVX__)
    c17_4 = _mm_add_pd(c17_4, _mm_mul_pd(a17_4, _mm256_castpd256_pd128(b17)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c17_4 = _mm_add_pd(c17_4, _mm_mul_pd(a17_4, b17));
#endif
    _mm_storeu_pd(&C[(l_n*56)+39], c17_4);
#else
    C[(l_n*56)+35] += A[168] * B[(l_n*21)+17];
    C[(l_n*56)+36] += A[169] * B[(l_n*21)+17];
    C[(l_n*56)+37] += A[170] * B[(l_n*21)+17];
    C[(l_n*56)+38] += A[171] * B[(l_n*21)+17];
    C[(l_n*56)+39] += A[172] * B[(l_n*21)+17];
    C[(l_n*56)+40] += A[173] * B[(l_n*21)+17];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b18 = _mm256_broadcast_sd(&B[(l_n*21)+18]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b18 = _mm_loaddup_pd(&B[(l_n*21)+18]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c18_0 = _mm256_loadu_pd(&C[(l_n*56)+35]);
    __m256d a18_0 = _mm256_loadu_pd(&A[174]);
    c18_0 = _mm256_add_pd(c18_0, _mm256_mul_pd(a18_0, b18));
    _mm256_storeu_pd(&C[(l_n*56)+35], c18_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c18_0 = _mm_loadu_pd(&C[(l_n*56)+35]);
    __m128d a18_0 = _mm_loadu_pd(&A[174]);
    c18_0 = _mm_add_pd(c18_0, _mm_mul_pd(a18_0, b18));
    _mm_storeu_pd(&C[(l_n*56)+35], c18_0);
    __m128d c18_2 = _mm_loadu_pd(&C[(l_n*56)+37]);
    __m128d a18_2 = _mm_loadu_pd(&A[176]);
    c18_2 = _mm_add_pd(c18_2, _mm_mul_pd(a18_2, b18));
    _mm_storeu_pd(&C[(l_n*56)+37], c18_2);
#endif
    __m128d c18_4 = _mm_loadu_pd(&C[(l_n*56)+39]);
    __m128d a18_4 = _mm_loadu_pd(&A[178]);
#if defined(__SSE3__) && defined(__AVX__)
    c18_4 = _mm_add_pd(c18_4, _mm_mul_pd(a18_4, _mm256_castpd256_pd128(b18)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c18_4 = _mm_add_pd(c18_4, _mm_mul_pd(a18_4, b18));
#endif
    _mm_storeu_pd(&C[(l_n*56)+39], c18_4);
#else
    C[(l_n*56)+35] += A[174] * B[(l_n*21)+18];
    C[(l_n*56)+36] += A[175] * B[(l_n*21)+18];
    C[(l_n*56)+37] += A[176] * B[(l_n*21)+18];
    C[(l_n*56)+38] += A[177] * B[(l_n*21)+18];
    C[(l_n*56)+39] += A[178] * B[(l_n*21)+18];
    C[(l_n*56)+40] += A[179] * B[(l_n*21)+18];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b19 = _mm256_broadcast_sd(&B[(l_n*21)+19]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b19 = _mm_loaddup_pd(&B[(l_n*21)+19]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c19_0 = _mm256_loadu_pd(&C[(l_n*56)+35]);
    __m256d a19_0 = _mm256_loadu_pd(&A[180]);
    c19_0 = _mm256_add_pd(c19_0, _mm256_mul_pd(a19_0, b19));
    _mm256_storeu_pd(&C[(l_n*56)+35], c19_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c19_0 = _mm_loadu_pd(&C[(l_n*56)+35]);
    __m128d a19_0 = _mm_loadu_pd(&A[180]);
    c19_0 = _mm_add_pd(c19_0, _mm_mul_pd(a19_0, b19));
    _mm_storeu_pd(&C[(l_n*56)+35], c19_0);
    __m128d c19_2 = _mm_loadu_pd(&C[(l_n*56)+37]);
    __m128d a19_2 = _mm_loadu_pd(&A[182]);
    c19_2 = _mm_add_pd(c19_2, _mm_mul_pd(a19_2, b19));
    _mm_storeu_pd(&C[(l_n*56)+37], c19_2);
#endif
    __m128d c19_4 = _mm_loadu_pd(&C[(l_n*56)+39]);
    __m128d a19_4 = _mm_loadu_pd(&A[184]);
#if defined(__SSE3__) && defined(__AVX__)
    c19_4 = _mm_add_pd(c19_4, _mm_mul_pd(a19_4, _mm256_castpd256_pd128(b19)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c19_4 = _mm_add_pd(c19_4, _mm_mul_pd(a19_4, b19));
#endif
    _mm_storeu_pd(&C[(l_n*56)+39], c19_4);
#else
    C[(l_n*56)+35] += A[180] * B[(l_n*21)+19];
    C[(l_n*56)+36] += A[181] * B[(l_n*21)+19];
    C[(l_n*56)+37] += A[182] * B[(l_n*21)+19];
    C[(l_n*56)+38] += A[183] * B[(l_n*21)+19];
    C[(l_n*56)+39] += A[184] * B[(l_n*21)+19];
    C[(l_n*56)+40] += A[185] * B[(l_n*21)+19];
#endif

#if defined(__SSE3__) || defined(__AVX__)
#if defined(__SSE3__) && defined(__AVX__)
    __m256d b20 = _mm256_broadcast_sd(&B[(l_n*21)+20]);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d b20 = _mm_loaddup_pd(&B[(l_n*21)+20]);
#endif
#if defined(__SSE3__) && defined(__AVX__)
    __m256d c20_0 = _mm256_loadu_pd(&C[(l_n*56)+35]);
    __m256d a20_0 = _mm256_loadu_pd(&A[186]);
    c20_0 = _mm256_add_pd(c20_0, _mm256_mul_pd(a20_0, b20));
    _mm256_storeu_pd(&C[(l_n*56)+35], c20_0);
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    __m128d c20_0 = _mm_loadu_pd(&C[(l_n*56)+35]);
    __m128d a20_0 = _mm_loadu_pd(&A[186]);
    c20_0 = _mm_add_pd(c20_0, _mm_mul_pd(a20_0, b20));
    _mm_storeu_pd(&C[(l_n*56)+35], c20_0);
    __m128d c20_2 = _mm_loadu_pd(&C[(l_n*56)+37]);
    __m128d a20_2 = _mm_loadu_pd(&A[188]);
    c20_2 = _mm_add_pd(c20_2, _mm_mul_pd(a20_2, b20));
    _mm_storeu_pd(&C[(l_n*56)+37], c20_2);
#endif
    __m128d c20_4 = _mm_loadu_pd(&C[(l_n*56)+39]);
    __m128d a20_4 = _mm_loadu_pd(&A[190]);
#if defined(__SSE3__) && defined(__AVX__)
    c20_4 = _mm_add_pd(c20_4, _mm_mul_pd(a20_4, _mm256_castpd256_pd128(b20)));
#endif
#if defined(__SSE3__) && !defined(__AVX__)
    c20_4 = _mm_add_pd(c20_4, _mm_mul_pd(a20_4, b20));
#endif
    _mm_storeu_pd(&C[(l_n*56)+39], c20_4);
#else
    C[(l_n*56)+35] += A[186] * B[(l_n*21)+20];
    C[(l_n*56)+36] += A[187] * B[(l_n*21)+20];
    C[(l_n*56)+37] += A[188] * B[(l_n*21)+20];
    C[(l_n*56)+38] += A[189] * B[(l_n*21)+20];
    C[(l_n*56)+39] += A[190] * B[(l_n*21)+20];
    C[(l_n*56)+40] += A[191] * B[(l_n*21)+20];
#endif

  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 3456;
#endif
}

void libxsmm_m56_n9_k21_ldA56_ldB21_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vmovapd 0(%%rdx), %%ymm4\n\t"
                       "vmovapd 32(%%rdx), %%ymm5\n\t"
                       "vmovapd 64(%%rdx), %%ymm6\n\t"
                       "vmovapd 96(%%rdx), %%ymm7\n\t"
                       "vmovapd 448(%%rdx), %%ymm8\n\t"
                       "vmovapd 480(%%rdx), %%ymm9\n\t"
                       "vmovapd 512(%%rdx), %%ymm10\n\t"
                       "vmovapd 544(%%rdx), %%ymm11\n\t"
                       "vmovapd 896(%%rdx), %%ymm12\n\t"
                       "vmovapd 928(%%rdx), %%ymm13\n\t"
                       "vmovapd 960(%%rdx), %%ymm14\n\t"
                       "vmovapd 992(%%rdx), %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 168(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 336(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 176(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 344(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 184(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 352(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 360(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 368(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 376(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 216(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 384(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 224(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 392(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 232(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 400(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 240(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 408(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 248(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 416(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 256(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 424(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 264(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 432(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 272(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 440(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 280(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 456(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 296(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 464(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 304(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 472(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 312(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 480(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 320(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 488(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 328(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 496(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm4, 0(%%rdx)\n\t"
                       "vmovapd %%ymm5, 32(%%rdx)\n\t"
                       "vmovapd %%ymm6, 64(%%rdx)\n\t"
                       "vmovapd %%ymm7, 96(%%rdx)\n\t"
                       "vmovapd %%ymm8, 448(%%rdx)\n\t"
                       "vmovapd %%ymm9, 480(%%rdx)\n\t"
                       "vmovapd %%ymm10, 512(%%rdx)\n\t"
                       "vmovapd %%ymm11, 544(%%rdx)\n\t"
                       "vmovapd %%ymm12, 896(%%rdx)\n\t"
                       "vmovapd %%ymm13, 928(%%rdx)\n\t"
                       "vmovapd %%ymm14, 960(%%rdx)\n\t"
                       "vmovapd %%ymm15, 992(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $9280, %%rdi\n\t"
                       "cmpq $32, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $32, %%r10\n\t"
                       "34:\n\t"
                       "addq $12, %%r10\n\t"
                       "vmovapd 0(%%rdx), %%ymm7\n\t"
                       "vmovapd 32(%%rdx), %%ymm8\n\t"
                       "vmovapd 64(%%rdx), %%ymm9\n\t"
                       "vmovapd 448(%%rdx), %%ymm10\n\t"
                       "vmovapd 480(%%rdx), %%ymm11\n\t"
                       "vmovapd 512(%%rdx), %%ymm12\n\t"
                       "vmovapd 896(%%rdx), %%ymm13\n\t"
                       "vmovapd 928(%%rdx), %%ymm14\n\t"
                       "vmovapd 960(%%rdx), %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 168(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 336(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 176(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 344(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 184(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 352(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 360(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 368(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 376(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 216(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 384(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 224(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 392(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 232(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 400(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 240(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 408(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 248(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 416(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 256(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 424(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 264(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 432(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 272(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 440(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 280(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 456(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 296(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 464(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 304(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 472(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 312(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 480(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 320(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 488(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 328(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 496(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm7, 0(%%rdx)\n\t"
                       "vmovapd %%ymm8, 32(%%rdx)\n\t"
                       "vmovapd %%ymm9, 64(%%rdx)\n\t"
                       "vmovapd %%ymm10, 448(%%rdx)\n\t"
                       "vmovapd %%ymm11, 480(%%rdx)\n\t"
                       "vmovapd %%ymm12, 512(%%rdx)\n\t"
                       "vmovapd %%ymm13, 896(%%rdx)\n\t"
                       "vmovapd %%ymm14, 928(%%rdx)\n\t"
                       "vmovapd %%ymm15, 960(%%rdx)\n\t"
                       "addq $96, %%rdx\n\t"
                       "subq $9312, %%rdi\n\t"
                       "cmpq $56, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $896, %%rdx\n\t"
                       "addq $504, %%rsi\n\t"
                       "subq $448, %%rdi\n\t"
                       "cmpq $9, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 21168;
#endif
}

void libxsmm_m36_n9_k55_ldA36_ldB56_ldC36_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vxorpd %%ymm4, %%ymm4, %%ymm4\n\t"
                       "vxorpd %%ymm5, %%ymm5, %%ymm5\n\t"
                       "vxorpd %%ymm6, %%ymm6, %%ymm6\n\t"
                       "vxorpd %%ymm7, %%ymm7, %%ymm7\n\t"
                       "vxorpd %%ymm8, %%ymm8, %%ymm8\n\t"
                       "vxorpd %%ymm9, %%ymm9, %%ymm9\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "cmpq $52, %%r12\n\t"
                       "jl 35b\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "subq $440, %%rsi\n\t"
                       "vmovapd %%ymm4, 0(%%rdx)\n\t"
                       "vmovapd %%ymm5, 32(%%rdx)\n\t"
                       "vmovapd %%ymm6, 64(%%rdx)\n\t"
                       "vmovapd %%ymm7, 96(%%rdx)\n\t"
                       "vmovapd %%ymm8, 288(%%rdx)\n\t"
                       "vmovapd %%ymm9, 320(%%rdx)\n\t"
                       "vmovapd %%ymm10, 352(%%rdx)\n\t"
                       "vmovapd %%ymm11, 384(%%rdx)\n\t"
                       "vmovapd %%ymm12, 576(%%rdx)\n\t"
                       "vmovapd %%ymm13, 608(%%rdx)\n\t"
                       "vmovapd %%ymm14, 640(%%rdx)\n\t"
                       "vmovapd %%ymm15, 672(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $15712, %%rdi\n\t"
                       "cmpq $32, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $32, %%r10\n\t"
                       "34:\n\t"
                       "addq $4, %%r10\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "cmpq $52, %%r12\n\t"
                       "jl 35b\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "subq $440, %%rsi\n\t"
                       "vmovapd %%ymm13, 0(%%rdx)\n\t"
                       "vmovapd %%ymm14, 288(%%rdx)\n\t"
                       "vmovapd %%ymm15, 576(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $15808, %%rdi\n\t"
                       "cmpq $36, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $576, %%rdx\n\t"
                       "addq $1344, %%rsi\n\t"
                       "subq $288, %%rdi\n\t"
                       "cmpq $9, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 35640;
#endif
}

void libxsmm_m36_n9_k9_ldA36_ldB9_ldC36_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vxorpd %%ymm4, %%ymm4, %%ymm4\n\t"
                       "vxorpd %%ymm5, %%ymm5, %%ymm5\n\t"
                       "vxorpd %%ymm6, %%ymm6, %%ymm6\n\t"
                       "vxorpd %%ymm7, %%ymm7, %%ymm7\n\t"
                       "vxorpd %%ymm8, %%ymm8, %%ymm8\n\t"
                       "vxorpd %%ymm9, %%ymm9, %%ymm9\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 168(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 176(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 184(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm4, 0(%%rdx)\n\t"
                       "vmovapd %%ymm5, 32(%%rdx)\n\t"
                       "vmovapd %%ymm6, 64(%%rdx)\n\t"
                       "vmovapd %%ymm7, 96(%%rdx)\n\t"
                       "vmovapd %%ymm8, 288(%%rdx)\n\t"
                       "vmovapd %%ymm9, 320(%%rdx)\n\t"
                       "vmovapd %%ymm10, 352(%%rdx)\n\t"
                       "vmovapd %%ymm11, 384(%%rdx)\n\t"
                       "vmovapd %%ymm12, 576(%%rdx)\n\t"
                       "vmovapd %%ymm13, 608(%%rdx)\n\t"
                       "vmovapd %%ymm14, 640(%%rdx)\n\t"
                       "vmovapd %%ymm15, 672(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $2464, %%rdi\n\t"
                       "cmpq $32, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $32, %%r10\n\t"
                       "34:\n\t"
                       "addq $4, %%r10\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 168(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 176(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 184(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm13, 0(%%rdx)\n\t"
                       "vmovapd %%ymm14, 288(%%rdx)\n\t"
                       "vmovapd %%ymm15, 576(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $2560, %%rdi\n\t"
                       "cmpq $36, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $576, %%rdx\n\t"
                       "addq $216, %%rsi\n\t"
                       "subq $288, %%rdi\n\t"
                       "cmpq $9, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 5832;
#endif
}

void libxsmm_m36_n9_k9_ldA36_ldB9_ldC36_alpha1_beta1_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vmovapd 0(%%rdx), %%ymm4\n\t"
                       "vmovapd 32(%%rdx), %%ymm5\n\t"
                       "vmovapd 64(%%rdx), %%ymm6\n\t"
                       "vmovapd 96(%%rdx), %%ymm7\n\t"
                       "vmovapd 288(%%rdx), %%ymm8\n\t"
                       "vmovapd 320(%%rdx), %%ymm9\n\t"
                       "vmovapd 352(%%rdx), %%ymm10\n\t"
                       "vmovapd 384(%%rdx), %%ymm11\n\t"
                       "vmovapd 576(%%rdx), %%ymm12\n\t"
                       "vmovapd 608(%%rdx), %%ymm13\n\t"
                       "vmovapd 640(%%rdx), %%ymm14\n\t"
                       "vmovapd 672(%%rdx), %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 168(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 176(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 184(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm4, 0(%%rdx)\n\t"
                       "vmovapd %%ymm5, 32(%%rdx)\n\t"
                       "vmovapd %%ymm6, 64(%%rdx)\n\t"
                       "vmovapd %%ymm7, 96(%%rdx)\n\t"
                       "vmovapd %%ymm8, 288(%%rdx)\n\t"
                       "vmovapd %%ymm9, 320(%%rdx)\n\t"
                       "vmovapd %%ymm10, 352(%%rdx)\n\t"
                       "vmovapd %%ymm11, 384(%%rdx)\n\t"
                       "vmovapd %%ymm12, 576(%%rdx)\n\t"
                       "vmovapd %%ymm13, 608(%%rdx)\n\t"
                       "vmovapd %%ymm14, 640(%%rdx)\n\t"
                       "vmovapd %%ymm15, 672(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $2464, %%rdi\n\t"
                       "cmpq $32, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $32, %%r10\n\t"
                       "34:\n\t"
                       "addq $4, %%r10\n\t"
                       "vmovapd 0(%%rdx), %%ymm13\n\t"
                       "vmovapd 288(%%rdx), %%ymm14\n\t"
                       "vmovapd 576(%%rdx), %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 168(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 176(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 184(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm13, 0(%%rdx)\n\t"
                       "vmovapd %%ymm14, 288(%%rdx)\n\t"
                       "vmovapd %%ymm15, 576(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $2560, %%rdi\n\t"
                       "cmpq $36, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $576, %%rdx\n\t"
                       "addq $216, %%rsi\n\t"
                       "subq $288, %%rdi\n\t"
                       "cmpq $9, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 5832;
#endif
}

void libxsmm_m20_n9_k34_ldA36_ldB36_ldC20_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vxorpd %%ymm4, %%ymm4, %%ymm4\n\t"
                       "vxorpd %%ymm5, %%ymm5, %%ymm5\n\t"
                       "vxorpd %%ymm6, %%ymm6, %%ymm6\n\t"
                       "vxorpd %%ymm7, %%ymm7, %%ymm7\n\t"
                       "vxorpd %%ymm8, %%ymm8, %%ymm8\n\t"
                       "vxorpd %%ymm9, %%ymm9, %%ymm9\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 576(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 576(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 576(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 576(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "cmpq $32, %%r12\n\t"
                       "jl 35b\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 576(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 576(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "subq $272, %%rsi\n\t"
                       "vmovapd %%ymm4, 0(%%rdx)\n\t"
                       "vmovapd %%ymm5, 32(%%rdx)\n\t"
                       "vmovapd %%ymm6, 64(%%rdx)\n\t"
                       "vmovapd %%ymm7, 96(%%rdx)\n\t"
                       "vmovapd %%ymm8, 160(%%rdx)\n\t"
                       "vmovapd %%ymm9, 192(%%rdx)\n\t"
                       "vmovapd %%ymm10, 224(%%rdx)\n\t"
                       "vmovapd %%ymm11, 256(%%rdx)\n\t"
                       "vmovapd %%ymm12, 320(%%rdx)\n\t"
                       "vmovapd %%ymm13, 352(%%rdx)\n\t"
                       "vmovapd %%ymm14, 384(%%rdx)\n\t"
                       "vmovapd %%ymm15, 416(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $9664, %%rdi\n\t"
                       "cmpq $16, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $16, %%r10\n\t"
                       "34:\n\t"
                       "addq $4, %%r10\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 576(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 576(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 576(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 576(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "cmpq $32, %%r12\n\t"
                       "jl 35b\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 576(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 576(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "subq $272, %%rsi\n\t"
                       "vmovapd %%ymm13, 0(%%rdx)\n\t"
                       "vmovapd %%ymm14, 160(%%rdx)\n\t"
                       "vmovapd %%ymm15, 320(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $9760, %%rdi\n\t"
                       "cmpq $20, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $320, %%rdx\n\t"
                       "addq $864, %%rsi\n\t"
                       "subq $160, %%rdi\n\t"
                       "cmpq $9, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 12240;
#endif
}

void libxsmm_m20_n9_k9_ldA20_ldB9_ldC20_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vxorpd %%ymm4, %%ymm4, %%ymm4\n\t"
                       "vxorpd %%ymm5, %%ymm5, %%ymm5\n\t"
                       "vxorpd %%ymm6, %%ymm6, %%ymm6\n\t"
                       "vxorpd %%ymm7, %%ymm7, %%ymm7\n\t"
                       "vxorpd %%ymm8, %%ymm8, %%ymm8\n\t"
                       "vxorpd %%ymm9, %%ymm9, %%ymm9\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 168(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 176(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 184(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm4, 0(%%rdx)\n\t"
                       "vmovapd %%ymm5, 32(%%rdx)\n\t"
                       "vmovapd %%ymm6, 64(%%rdx)\n\t"
                       "vmovapd %%ymm7, 96(%%rdx)\n\t"
                       "vmovapd %%ymm8, 160(%%rdx)\n\t"
                       "vmovapd %%ymm9, 192(%%rdx)\n\t"
                       "vmovapd %%ymm10, 224(%%rdx)\n\t"
                       "vmovapd %%ymm11, 256(%%rdx)\n\t"
                       "vmovapd %%ymm12, 320(%%rdx)\n\t"
                       "vmovapd %%ymm13, 352(%%rdx)\n\t"
                       "vmovapd %%ymm14, 384(%%rdx)\n\t"
                       "vmovapd %%ymm15, 416(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $1312, %%rdi\n\t"
                       "cmpq $16, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $16, %%r10\n\t"
                       "34:\n\t"
                       "addq $4, %%r10\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 168(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 176(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 184(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm13, 0(%%rdx)\n\t"
                       "vmovapd %%ymm14, 160(%%rdx)\n\t"
                       "vmovapd %%ymm15, 320(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $1408, %%rdi\n\t"
                       "cmpq $20, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $320, %%rdx\n\t"
                       "addq $216, %%rsi\n\t"
                       "subq $160, %%rdi\n\t"
                       "cmpq $9, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 3240;
#endif
}

void libxsmm_m20_n9_k9_ldA20_ldB9_ldC20_alpha1_beta1_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vmovapd 0(%%rdx), %%ymm4\n\t"
                       "vmovapd 32(%%rdx), %%ymm5\n\t"
                       "vmovapd 64(%%rdx), %%ymm6\n\t"
                       "vmovapd 96(%%rdx), %%ymm7\n\t"
                       "vmovapd 160(%%rdx), %%ymm8\n\t"
                       "vmovapd 192(%%rdx), %%ymm9\n\t"
                       "vmovapd 224(%%rdx), %%ymm10\n\t"
                       "vmovapd 256(%%rdx), %%ymm11\n\t"
                       "vmovapd 320(%%rdx), %%ymm12\n\t"
                       "vmovapd 352(%%rdx), %%ymm13\n\t"
                       "vmovapd 384(%%rdx), %%ymm14\n\t"
                       "vmovapd 416(%%rdx), %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 168(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 176(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 184(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm4, 0(%%rdx)\n\t"
                       "vmovapd %%ymm5, 32(%%rdx)\n\t"
                       "vmovapd %%ymm6, 64(%%rdx)\n\t"
                       "vmovapd %%ymm7, 96(%%rdx)\n\t"
                       "vmovapd %%ymm8, 160(%%rdx)\n\t"
                       "vmovapd %%ymm9, 192(%%rdx)\n\t"
                       "vmovapd %%ymm10, 224(%%rdx)\n\t"
                       "vmovapd %%ymm11, 256(%%rdx)\n\t"
                       "vmovapd %%ymm12, 320(%%rdx)\n\t"
                       "vmovapd %%ymm13, 352(%%rdx)\n\t"
                       "vmovapd %%ymm14, 384(%%rdx)\n\t"
                       "vmovapd %%ymm15, 416(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $1312, %%rdi\n\t"
                       "cmpq $16, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $16, %%r10\n\t"
                       "34:\n\t"
                       "addq $4, %%r10\n\t"
                       "vmovapd 0(%%rdx), %%ymm13\n\t"
                       "vmovapd 160(%%rdx), %%ymm14\n\t"
                       "vmovapd 320(%%rdx), %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 168(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 176(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 184(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $160, %%rdi\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm13, 0(%%rdx)\n\t"
                       "vmovapd %%ymm14, 160(%%rdx)\n\t"
                       "vmovapd %%ymm15, 320(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $1408, %%rdi\n\t"
                       "cmpq $20, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $320, %%rdx\n\t"
                       "addq $216, %%rsi\n\t"
                       "subq $160, %%rdi\n\t"
                       "cmpq $9, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 3240;
#endif
}

void libxsmm_m12_n9_k19_ldA36_ldB20_ldC12_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $12, %%r10\n\t"
                       "vxorpd %%ymm7, %%ymm7, %%ymm7\n\t"
                       "vxorpd %%ymm8, %%ymm8, %%ymm8\n\t"
                       "vxorpd %%ymm9, %%ymm9, %%ymm9\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 320(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 168(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 328(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 176(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 336(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 184(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 344(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 352(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 360(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 368(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 216(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 376(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 224(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 384(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 232(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 392(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 240(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 400(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 248(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 408(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 256(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 416(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 264(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 424(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 272(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 432(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 280(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 440(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 296(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 456(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 304(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 464(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $288, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm7, 0(%%rdx)\n\t"
                       "vmovapd %%ymm8, 32(%%rdx)\n\t"
                       "vmovapd %%ymm9, 64(%%rdx)\n\t"
                       "vmovapd %%ymm10, 96(%%rdx)\n\t"
                       "vmovapd %%ymm11, 128(%%rdx)\n\t"
                       "vmovapd %%ymm12, 160(%%rdx)\n\t"
                       "vmovapd %%ymm13, 192(%%rdx)\n\t"
                       "vmovapd %%ymm14, 224(%%rdx)\n\t"
                       "vmovapd %%ymm15, 256(%%rdx)\n\t"
                       "addq $96, %%rdx\n\t"
                       "subq $5376, %%rdi\n\t"
                       "cmpq $12, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $192, %%rdx\n\t"
                       "addq $480, %%rsi\n\t"
                       "subq $96, %%rdi\n\t"
                       "cmpq $9, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 4104;
#endif
}

void libxsmm_m12_n9_k9_ldA12_ldB9_ldC12_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $12, %%r10\n\t"
                       "vxorpd %%ymm7, %%ymm7, %%ymm7\n\t"
                       "vxorpd %%ymm8, %%ymm8, %%ymm8\n\t"
                       "vxorpd %%ymm9, %%ymm9, %%ymm9\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $96, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $96, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $96, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 168(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $96, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 176(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $96, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 184(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $96, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $96, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $96, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $96, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm7, 0(%%rdx)\n\t"
                       "vmovapd %%ymm8, 32(%%rdx)\n\t"
                       "vmovapd %%ymm9, 64(%%rdx)\n\t"
                       "vmovapd %%ymm10, 96(%%rdx)\n\t"
                       "vmovapd %%ymm11, 128(%%rdx)\n\t"
                       "vmovapd %%ymm12, 160(%%rdx)\n\t"
                       "vmovapd %%ymm13, 192(%%rdx)\n\t"
                       "vmovapd %%ymm14, 224(%%rdx)\n\t"
                       "vmovapd %%ymm15, 256(%%rdx)\n\t"
                       "addq $96, %%rdx\n\t"
                       "subq $768, %%rdi\n\t"
                       "cmpq $12, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $192, %%rdx\n\t"
                       "addq $216, %%rsi\n\t"
                       "subq $96, %%rdi\n\t"
                       "cmpq $9, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1944;
#endif
}

void libxsmm_m12_n9_k9_ldA12_ldB9_ldC12_alpha1_beta1_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $12, %%r10\n\t"
                       "vmovapd 0(%%rdx), %%ymm7\n\t"
                       "vmovapd 32(%%rdx), %%ymm8\n\t"
                       "vmovapd 64(%%rdx), %%ymm9\n\t"
                       "vmovapd 96(%%rdx), %%ymm10\n\t"
                       "vmovapd 128(%%rdx), %%ymm11\n\t"
                       "vmovapd 160(%%rdx), %%ymm12\n\t"
                       "vmovapd 192(%%rdx), %%ymm13\n\t"
                       "vmovapd 224(%%rdx), %%ymm14\n\t"
                       "vmovapd 256(%%rdx), %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $96, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $96, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $96, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 168(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $96, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 176(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $96, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 184(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $96, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $96, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $96, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $96, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm7, 0(%%rdx)\n\t"
                       "vmovapd %%ymm8, 32(%%rdx)\n\t"
                       "vmovapd %%ymm9, 64(%%rdx)\n\t"
                       "vmovapd %%ymm10, 96(%%rdx)\n\t"
                       "vmovapd %%ymm11, 128(%%rdx)\n\t"
                       "vmovapd %%ymm12, 160(%%rdx)\n\t"
                       "vmovapd %%ymm13, 192(%%rdx)\n\t"
                       "vmovapd %%ymm14, 224(%%rdx)\n\t"
                       "vmovapd %%ymm15, 256(%%rdx)\n\t"
                       "addq $96, %%rdx\n\t"
                       "subq $768, %%rdi\n\t"
                       "cmpq $12, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $192, %%rdx\n\t"
                       "addq $216, %%rsi\n\t"
                       "subq $96, %%rdi\n\t"
                       "cmpq $9, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1944;
#endif
}

void libxsmm_m4_n9_k9_ldA36_ldB12_ldC4_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $4, %%r10\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 216(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 224(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 232(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 240(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 248(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 256(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm13, 0(%%rdx)\n\t"
                       "vmovapd %%ymm14, 32(%%rdx)\n\t"
                       "vmovapd %%ymm15, 64(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $2560, %%rdi\n\t"
                       "cmpq $4, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $64, %%rdx\n\t"
                       "addq $288, %%rsi\n\t"
                       "subq $32, %%rdi\n\t"
                       "cmpq $9, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 648;
#endif
}

void libxsmm_m4_n9_k9_ldA4_ldB9_ldC4_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $4, %%r10\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 168(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 176(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 184(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm13, 0(%%rdx)\n\t"
                       "vmovapd %%ymm14, 32(%%rdx)\n\t"
                       "vmovapd %%ymm15, 64(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $256, %%rdi\n\t"
                       "cmpq $4, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $64, %%rdx\n\t"
                       "addq $216, %%rsi\n\t"
                       "subq $32, %%rdi\n\t"
                       "cmpq $9, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 648;
#endif
}

void libxsmm_m4_n9_k9_ldA4_ldB9_ldC4_alpha1_beta1_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $4, %%r10\n\t"
                       "vmovapd 0(%%rdx), %%ymm13\n\t"
                       "vmovapd 32(%%rdx), %%ymm14\n\t"
                       "vmovapd 64(%%rdx), %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 168(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 176(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 184(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm13, 0(%%rdx)\n\t"
                       "vmovapd %%ymm14, 32(%%rdx)\n\t"
                       "vmovapd %%ymm15, 64(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $256, %%rdi\n\t"
                       "cmpq $4, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $64, %%rdx\n\t"
                       "addq $216, %%rsi\n\t"
                       "subq $32, %%rdi\n\t"
                       "cmpq $9, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 648;
#endif
}

void libxsmm_m4_n9_k3_ldA36_ldB4_ldC4_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $4, %%r10\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $288, %%rdi\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm13, 0(%%rdx)\n\t"
                       "vmovapd %%ymm14, 32(%%rdx)\n\t"
                       "vmovapd %%ymm15, 64(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $832, %%rdi\n\t"
                       "cmpq $4, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $64, %%rdx\n\t"
                       "addq $96, %%rsi\n\t"
                       "subq $32, %%rdi\n\t"
                       "cmpq $9, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 216;
#endif
}

void libxsmm_m4_n21_k21_ldA4_ldB24_ldC4_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $4, %%r10\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 384(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 392(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 400(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 216(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 408(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 224(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 416(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 232(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 424(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 240(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 432(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 248(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 440(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 256(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 264(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 456(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 272(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 464(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 280(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 472(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 480(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 296(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 488(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 304(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 496(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 312(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 504(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 320(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 512(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 328(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 520(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 336(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 528(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 344(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 536(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 352(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 544(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm13, 0(%%rdx)\n\t"
                       "vmovapd %%ymm14, 32(%%rdx)\n\t"
                       "vmovapd %%ymm15, 64(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $640, %%rdi\n\t"
                       "cmpq $4, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $64, %%rdx\n\t"
                       "addq $576, %%rsi\n\t"
                       "subq $32, %%rdi\n\t"
                       "cmpq $21, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 3528;
#endif
}

void libxsmm_m4_n3_k21_ldA4_ldB24_ldC4_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $4, %%r10\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 384(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 392(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 400(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 216(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 408(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 224(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 416(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 232(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 424(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 240(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 432(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 248(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 440(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 256(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 264(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 456(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 272(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 464(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 280(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 472(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 480(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 296(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 488(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 304(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 496(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 312(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 504(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 320(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 512(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 328(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 520(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 336(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 528(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 344(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 536(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 352(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 544(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm13, 0(%%rdx)\n\t"
                       "vmovapd %%ymm14, 32(%%rdx)\n\t"
                       "vmovapd %%ymm15, 64(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $640, %%rdi\n\t"
                       "cmpq $4, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $64, %%rdx\n\t"
                       "addq $576, %%rsi\n\t"
                       "subq $32, %%rdi\n\t"
                       "cmpq $3, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 504;
#endif
}

void libxsmm_m24_n3_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vxorpd %%ymm4, %%ymm4, %%ymm4\n\t"
                       "vxorpd %%ymm5, %%ymm5, %%ymm5\n\t"
                       "vxorpd %%ymm6, %%ymm6, %%ymm6\n\t"
                       "vxorpd %%ymm7, %%ymm7, %%ymm7\n\t"
                       "vxorpd %%ymm8, %%ymm8, %%ymm8\n\t"
                       "vxorpd %%ymm9, %%ymm9, %%ymm9\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 384(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 392(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 400(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 216(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 408(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 224(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 416(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 232(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 424(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 240(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 432(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 248(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 440(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 256(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 264(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 456(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 272(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 464(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 280(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 472(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 480(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 296(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 488(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 304(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 496(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 312(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 504(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 320(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 512(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 328(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 520(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 336(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 528(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 344(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 536(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 352(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 544(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm4, 0(%%rdx)\n\t"
                       "vmovapd %%ymm5, 32(%%rdx)\n\t"
                       "vmovapd %%ymm6, 64(%%rdx)\n\t"
                       "vmovapd %%ymm7, 96(%%rdx)\n\t"
                       "vmovapd %%ymm8, 192(%%rdx)\n\t"
                       "vmovapd %%ymm9, 224(%%rdx)\n\t"
                       "vmovapd %%ymm10, 256(%%rdx)\n\t"
                       "vmovapd %%ymm11, 288(%%rdx)\n\t"
                       "vmovapd %%ymm12, 384(%%rdx)\n\t"
                       "vmovapd %%ymm13, 416(%%rdx)\n\t"
                       "vmovapd %%ymm14, 448(%%rdx)\n\t"
                       "vmovapd %%ymm15, 480(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $3904, %%rdi\n\t"
                       "cmpq $16, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $16, %%r10\n\t"
                       "34:\n\t"
                       "addq $8, %%r10\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 384(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 392(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 400(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 216(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 408(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 224(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 416(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 232(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 424(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 240(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 432(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 248(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 440(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 256(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 264(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 456(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 272(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 464(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 280(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 472(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 480(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 296(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 488(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 304(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 496(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 312(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 504(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 320(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 512(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 328(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 520(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 336(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 528(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 344(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 536(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 352(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 544(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm10, 0(%%rdx)\n\t"
                       "vmovapd %%ymm11, 32(%%rdx)\n\t"
                       "vmovapd %%ymm12, 192(%%rdx)\n\t"
                       "vmovapd %%ymm13, 224(%%rdx)\n\t"
                       "vmovapd %%ymm14, 384(%%rdx)\n\t"
                       "vmovapd %%ymm15, 416(%%rdx)\n\t"
                       "addq $64, %%rdx\n\t"
                       "subq $3968, %%rdi\n\t"
                       "cmpq $24, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $384, %%rdx\n\t"
                       "addq $576, %%rsi\n\t"
                       "subq $192, %%rdi\n\t"
                       "cmpq $3, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 3024;
#endif
}

void libxsmm_m16_n3_k21_ldA16_ldB24_ldC16_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vxorpd %%ymm4, %%ymm4, %%ymm4\n\t"
                       "vxorpd %%ymm5, %%ymm5, %%ymm5\n\t"
                       "vxorpd %%ymm6, %%ymm6, %%ymm6\n\t"
                       "vxorpd %%ymm7, %%ymm7, %%ymm7\n\t"
                       "vxorpd %%ymm8, %%ymm8, %%ymm8\n\t"
                       "vxorpd %%ymm9, %%ymm9, %%ymm9\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 384(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 392(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 400(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 216(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 408(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 224(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 416(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 232(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 424(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 240(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 432(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 248(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 440(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 256(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 264(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 456(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 272(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 464(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 280(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 472(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 480(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 296(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 488(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 304(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 496(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 312(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 504(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 320(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 512(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 328(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 520(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 336(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 528(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 344(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 536(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 352(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 544(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm4, 0(%%rdx)\n\t"
                       "vmovapd %%ymm5, 32(%%rdx)\n\t"
                       "vmovapd %%ymm6, 64(%%rdx)\n\t"
                       "vmovapd %%ymm7, 96(%%rdx)\n\t"
                       "vmovapd %%ymm8, 128(%%rdx)\n\t"
                       "vmovapd %%ymm9, 160(%%rdx)\n\t"
                       "vmovapd %%ymm10, 192(%%rdx)\n\t"
                       "vmovapd %%ymm11, 224(%%rdx)\n\t"
                       "vmovapd %%ymm12, 256(%%rdx)\n\t"
                       "vmovapd %%ymm13, 288(%%rdx)\n\t"
                       "vmovapd %%ymm14, 320(%%rdx)\n\t"
                       "vmovapd %%ymm15, 352(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $2560, %%rdi\n\t"
                       "cmpq $16, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $256, %%rdx\n\t"
                       "addq $576, %%rsi\n\t"
                       "subq $128, %%rdi\n\t"
                       "cmpq $3, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 2016;
#endif
}

void libxsmm_m64_n3_k21_ldA64_ldB24_ldC64_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vxorpd %%ymm4, %%ymm4, %%ymm4\n\t"
                       "vxorpd %%ymm5, %%ymm5, %%ymm5\n\t"
                       "vxorpd %%ymm6, %%ymm6, %%ymm6\n\t"
                       "vxorpd %%ymm7, %%ymm7, %%ymm7\n\t"
                       "vxorpd %%ymm8, %%ymm8, %%ymm8\n\t"
                       "vxorpd %%ymm9, %%ymm9, %%ymm9\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 384(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $512, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 392(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $512, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 400(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $512, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 216(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 408(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $512, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 224(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 416(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $512, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 232(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 424(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $512, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 240(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 432(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $512, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 248(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 440(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $512, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 256(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $512, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 264(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 456(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $512, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 272(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 464(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $512, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 280(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 472(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $512, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 288(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 480(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $512, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 296(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 488(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $512, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 304(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 496(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $512, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 312(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 504(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $512, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 320(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 512(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $512, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 328(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 520(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $512, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 336(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 528(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $512, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 344(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 536(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $512, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 352(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 544(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $512, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm4, 0(%%rdx)\n\t"
                       "vmovapd %%ymm5, 32(%%rdx)\n\t"
                       "vmovapd %%ymm6, 64(%%rdx)\n\t"
                       "vmovapd %%ymm7, 96(%%rdx)\n\t"
                       "vmovapd %%ymm8, 512(%%rdx)\n\t"
                       "vmovapd %%ymm9, 544(%%rdx)\n\t"
                       "vmovapd %%ymm10, 576(%%rdx)\n\t"
                       "vmovapd %%ymm11, 608(%%rdx)\n\t"
                       "vmovapd %%ymm12, 1024(%%rdx)\n\t"
                       "vmovapd %%ymm13, 1056(%%rdx)\n\t"
                       "vmovapd %%ymm14, 1088(%%rdx)\n\t"
                       "vmovapd %%ymm15, 1120(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $10624, %%rdi\n\t"
                       "cmpq $64, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $1024, %%rdx\n\t"
                       "addq $576, %%rsi\n\t"
                       "subq $512, %%rdi\n\t"
                       "cmpq $3, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 8064;
#endif
}

void libxsmm_m4_n3_k56_ldA4_ldB56_ldC4_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $4, %%r10\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $32, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "cmpq $56, %%r12\n\t"
                       "jl 35b\n\t"
                       "subq $448, %%rsi\n\t"
                       "vmovapd %%ymm13, 0(%%rdx)\n\t"
                       "vmovapd %%ymm14, 32(%%rdx)\n\t"
                       "vmovapd %%ymm15, 64(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $1760, %%rdi\n\t"
                       "cmpq $4, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $64, %%rdx\n\t"
                       "addq $1344, %%rsi\n\t"
                       "subq $32, %%rdi\n\t"
                       "cmpq $3, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 1344;
#endif
}

void libxsmmsparse_e60658f2bd24617a41c922e81bd775a8_m4_n3_k3_ldA4_ldB0_ldC4_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_m = 0;
  unsigned int l_n = 0;
  for ( l_n = 0; l_n < 3; l_n++) {
    #pragma simd
    #pragma vector aligned
    for ( l_m = 0; l_m < 4; l_m++) { C[(l_n*4)+l_m] = 0.0; }
  }

  #pragma simd vectorlength(4)
  #pragma vector aligned
  for ( l_m = 0; l_m < 4; l_m++) {
    C[0+l_m] += A[0+l_m] * B[0];
    C[4+l_m] += A[4+l_m] * B[1];
    C[8+l_m] += A[8+l_m] * B[2];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 24;
#endif
}

void libxsmm_m16_n3_k56_ldA16_ldB56_ldC16_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vxorpd %%ymm4, %%ymm4, %%ymm4\n\t"
                       "vxorpd %%ymm5, %%ymm5, %%ymm5\n\t"
                       "vxorpd %%ymm6, %%ymm6, %%ymm6\n\t"
                       "vxorpd %%ymm7, %%ymm7, %%ymm7\n\t"
                       "vxorpd %%ymm8, %%ymm8, %%ymm8\n\t"
                       "vxorpd %%ymm9, %%ymm9, %%ymm9\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $128, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "cmpq $56, %%r12\n\t"
                       "jl 35b\n\t"
                       "subq $448, %%rsi\n\t"
                       "vmovapd %%ymm4, 0(%%rdx)\n\t"
                       "vmovapd %%ymm5, 32(%%rdx)\n\t"
                       "vmovapd %%ymm6, 64(%%rdx)\n\t"
                       "vmovapd %%ymm7, 96(%%rdx)\n\t"
                       "vmovapd %%ymm8, 128(%%rdx)\n\t"
                       "vmovapd %%ymm9, 160(%%rdx)\n\t"
                       "vmovapd %%ymm10, 192(%%rdx)\n\t"
                       "vmovapd %%ymm11, 224(%%rdx)\n\t"
                       "vmovapd %%ymm12, 256(%%rdx)\n\t"
                       "vmovapd %%ymm13, 288(%%rdx)\n\t"
                       "vmovapd %%ymm14, 320(%%rdx)\n\t"
                       "vmovapd %%ymm15, 352(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $7040, %%rdi\n\t"
                       "cmpq $16, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $256, %%rdx\n\t"
                       "addq $1344, %%rsi\n\t"
                       "subq $128, %%rdi\n\t"
                       "cmpq $3, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 5376;
#endif
}

void libxsmmsparse_e60658f2bd24617a41c922e81bd775a8_m16_n3_k3_ldA16_ldB0_ldC16_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_m = 0;
  unsigned int l_n = 0;
  for ( l_n = 0; l_n < 3; l_n++) {
    #pragma simd
    #pragma vector aligned
    for ( l_m = 0; l_m < 16; l_m++) { C[(l_n*16)+l_m] = 0.0; }
  }

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 16; l_m++) {
    C[0+l_m] += A[0+l_m] * B[0];
    C[16+l_m] += A[16+l_m] * B[1];
    C[32+l_m] += A[32+l_m] * B[2];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 96;
#endif
}

void libxsmmsparse_e60658f2bd24617a41c922e81bd775a8_m56_n3_k3_ldA56_ldB0_ldC56_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_m = 0;
  unsigned int l_n = 0;
  for ( l_n = 0; l_n < 3; l_n++) {
    #pragma simd
    #pragma vector aligned
    for ( l_m = 0; l_m < 56; l_m++) { C[(l_n*56)+l_m] = 0.0; }
  }

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 56; l_m++) {
    C[0+l_m] += A[0+l_m] * B[0];
    C[56+l_m] += A[56+l_m] * B[1];
    C[112+l_m] += A[112+l_m] * B[2];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 336;
#endif
}

void libxsmm_m64_n3_k56_ldA64_ldB56_ldC64_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vxorpd %%ymm4, %%ymm4, %%ymm4\n\t"
                       "vxorpd %%ymm5, %%ymm5, %%ymm5\n\t"
                       "vxorpd %%ymm6, %%ymm6, %%ymm6\n\t"
                       "vxorpd %%ymm7, %%ymm7, %%ymm7\n\t"
                       "vxorpd %%ymm8, %%ymm8, %%ymm8\n\t"
                       "vxorpd %%ymm9, %%ymm9, %%ymm9\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $512, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $512, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $512, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $512, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "cmpq $56, %%r12\n\t"
                       "jl 35b\n\t"
                       "subq $448, %%rsi\n\t"
                       "vmovapd %%ymm4, 0(%%rdx)\n\t"
                       "vmovapd %%ymm5, 32(%%rdx)\n\t"
                       "vmovapd %%ymm6, 64(%%rdx)\n\t"
                       "vmovapd %%ymm7, 96(%%rdx)\n\t"
                       "vmovapd %%ymm8, 512(%%rdx)\n\t"
                       "vmovapd %%ymm9, 544(%%rdx)\n\t"
                       "vmovapd %%ymm10, 576(%%rdx)\n\t"
                       "vmovapd %%ymm11, 608(%%rdx)\n\t"
                       "vmovapd %%ymm12, 1024(%%rdx)\n\t"
                       "vmovapd %%ymm13, 1056(%%rdx)\n\t"
                       "vmovapd %%ymm14, 1088(%%rdx)\n\t"
                       "vmovapd %%ymm15, 1120(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $28544, %%rdi\n\t"
                       "cmpq $64, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $1024, %%rdx\n\t"
                       "addq $1344, %%rsi\n\t"
                       "subq $512, %%rdi\n\t"
                       "cmpq $3, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 21504;
#endif
}

void libxsmm_m24_n3_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vxorpd %%ymm4, %%ymm4, %%ymm4\n\t"
                       "vxorpd %%ymm5, %%ymm5, %%ymm5\n\t"
                       "vxorpd %%ymm6, %%ymm6, %%ymm6\n\t"
                       "vxorpd %%ymm7, %%ymm7, %%ymm7\n\t"
                       "vxorpd %%ymm8, %%ymm8, %%ymm8\n\t"
                       "vxorpd %%ymm9, %%ymm9, %%ymm9\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "cmpq $56, %%r12\n\t"
                       "jl 35b\n\t"
                       "subq $448, %%rsi\n\t"
                       "vmovapd %%ymm4, 0(%%rdx)\n\t"
                       "vmovapd %%ymm5, 32(%%rdx)\n\t"
                       "vmovapd %%ymm6, 64(%%rdx)\n\t"
                       "vmovapd %%ymm7, 96(%%rdx)\n\t"
                       "vmovapd %%ymm8, 192(%%rdx)\n\t"
                       "vmovapd %%ymm9, 224(%%rdx)\n\t"
                       "vmovapd %%ymm10, 256(%%rdx)\n\t"
                       "vmovapd %%ymm11, 288(%%rdx)\n\t"
                       "vmovapd %%ymm12, 384(%%rdx)\n\t"
                       "vmovapd %%ymm13, 416(%%rdx)\n\t"
                       "vmovapd %%ymm14, 448(%%rdx)\n\t"
                       "vmovapd %%ymm15, 480(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $10624, %%rdi\n\t"
                       "cmpq $16, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $16, %%r10\n\t"
                       "34:\n\t"
                       "addq $8, %%r10\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "addq $192, %%rdi\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm15\n\t"
                       "cmpq $56, %%r12\n\t"
                       "jl 35b\n\t"
                       "subq $448, %%rsi\n\t"
                       "vmovapd %%ymm10, 0(%%rdx)\n\t"
                       "vmovapd %%ymm11, 32(%%rdx)\n\t"
                       "vmovapd %%ymm12, 192(%%rdx)\n\t"
                       "vmovapd %%ymm13, 224(%%rdx)\n\t"
                       "vmovapd %%ymm14, 384(%%rdx)\n\t"
                       "vmovapd %%ymm15, 416(%%rdx)\n\t"
                       "addq $64, %%rdx\n\t"
                       "subq $10688, %%rdi\n\t"
                       "cmpq $24, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $384, %%rdx\n\t"
                       "addq $1344, %%rsi\n\t"
                       "subq $192, %%rdi\n\t"
                       "cmpq $3, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 8064;
#endif
}

void libxsmmsparse_e60658f2bd24617a41c922e81bd775a8_m24_n3_k3_ldA24_ldB0_ldC24_alpha1_beta1_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
  unsigned int l_m = 0;

  #pragma simd vectorlength(8)
  #pragma vector aligned
  for ( l_m = 0; l_m < 24; l_m++) {
    C[0+l_m] += A[0+l_m] * B[0];
    C[24+l_m] += A[24+l_m] * B[1];
    C[48+l_m] += A[48+l_m] * B[2];
  }

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 144;
#endif
}

void libxsmm_m52_n9_k56_ldA52_ldB56_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vxorpd %%ymm4, %%ymm4, %%ymm4\n\t"
                       "vxorpd %%ymm5, %%ymm5, %%ymm5\n\t"
                       "vxorpd %%ymm6, %%ymm6, %%ymm6\n\t"
                       "vxorpd %%ymm7, %%ymm7, %%ymm7\n\t"
                       "vxorpd %%ymm8, %%ymm8, %%ymm8\n\t"
                       "vxorpd %%ymm9, %%ymm9, %%ymm9\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $416, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $416, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $416, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $416, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "cmpq $56, %%r12\n\t"
                       "jl 35b\n\t"
                       "subq $448, %%rsi\n\t"
                       "vmovapd %%ymm4, 0(%%rdx)\n\t"
                       "vmovapd %%ymm5, 32(%%rdx)\n\t"
                       "vmovapd %%ymm6, 64(%%rdx)\n\t"
                       "vmovapd %%ymm7, 96(%%rdx)\n\t"
                       "vmovapd %%ymm8, 416(%%rdx)\n\t"
                       "vmovapd %%ymm9, 448(%%rdx)\n\t"
                       "vmovapd %%ymm10, 480(%%rdx)\n\t"
                       "vmovapd %%ymm11, 512(%%rdx)\n\t"
                       "vmovapd %%ymm12, 832(%%rdx)\n\t"
                       "vmovapd %%ymm13, 864(%%rdx)\n\t"
                       "vmovapd %%ymm14, 896(%%rdx)\n\t"
                       "vmovapd %%ymm15, 928(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $23168, %%rdi\n\t"
                       "cmpq $48, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $48, %%r10\n\t"
                       "34:\n\t"
                       "addq $4, %%r10\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $416, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $416, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $416, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $416, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 448(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 896(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "cmpq $56, %%r12\n\t"
                       "jl 35b\n\t"
                       "subq $448, %%rsi\n\t"
                       "vmovapd %%ymm13, 0(%%rdx)\n\t"
                       "vmovapd %%ymm14, 416(%%rdx)\n\t"
                       "vmovapd %%ymm15, 832(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $23264, %%rdi\n\t"
                       "cmpq $52, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $832, %%rdx\n\t"
                       "addq $1344, %%rsi\n\t"
                       "subq $416, %%rdi\n\t"
                       "cmpq $9, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 52416;
#endif
}

void libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vxorpd %%ymm4, %%ymm4, %%ymm4\n\t"
                       "vxorpd %%ymm5, %%ymm5, %%ymm5\n\t"
                       "vxorpd %%ymm6, %%ymm6, %%ymm6\n\t"
                       "vxorpd %%ymm7, %%ymm7, %%ymm7\n\t"
                       "vxorpd %%ymm8, %%ymm8, %%ymm8\n\t"
                       "vxorpd %%ymm9, %%ymm9, %%ymm9\n\t"
                       "vxorpd %%ymm10, %%ymm10, %%ymm10\n\t"
                       "vxorpd %%ymm11, %%ymm11, %%ymm11\n\t"
                       "vxorpd %%ymm12, %%ymm12, %%ymm12\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $416, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $416, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $416, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 168(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $416, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 176(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $416, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 184(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $416, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $416, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $416, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm2\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $416, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm4, 0(%%rdx)\n\t"
                       "vmovapd %%ymm5, 32(%%rdx)\n\t"
                       "vmovapd %%ymm6, 64(%%rdx)\n\t"
                       "vmovapd %%ymm7, 96(%%rdx)\n\t"
                       "vmovapd %%ymm8, 416(%%rdx)\n\t"
                       "vmovapd %%ymm9, 448(%%rdx)\n\t"
                       "vmovapd %%ymm10, 480(%%rdx)\n\t"
                       "vmovapd %%ymm11, 512(%%rdx)\n\t"
                       "vmovapd %%ymm12, 832(%%rdx)\n\t"
                       "vmovapd %%ymm13, 864(%%rdx)\n\t"
                       "vmovapd %%ymm14, 896(%%rdx)\n\t"
                       "vmovapd %%ymm15, 928(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $3616, %%rdi\n\t"
                       "cmpq $48, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $48, %%r10\n\t"
                       "34:\n\t"
                       "addq $4, %%r10\n\t"
                       "vxorpd %%ymm13, %%ymm13, %%ymm13\n\t"
                       "vxorpd %%ymm14, %%ymm14, %%ymm14\n\t"
                       "vxorpd %%ymm15, %%ymm15, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $416, %%rdi\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 72(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 144(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $416, %%rdi\n\t"
                       "vbroadcastsd 8(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 80(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 152(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $416, %%rdi\n\t"
                       "vbroadcastsd 16(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 88(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 160(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $416, %%rdi\n\t"
                       "vbroadcastsd 24(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 96(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 168(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $416, %%rdi\n\t"
                       "vbroadcastsd 32(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 104(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 176(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $416, %%rdi\n\t"
                       "vbroadcastsd 40(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 112(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 184(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $416, %%rdi\n\t"
                       "vbroadcastsd 48(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 120(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 192(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $416, %%rdi\n\t"
                       "vbroadcastsd 56(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 128(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 200(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "addq $416, %%rdi\n\t"
                       "vbroadcastsd 64(%%rsi), %%ymm0\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm13\n\t"
                       "vbroadcastsd 136(%%rsi), %%ymm1\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm14\n\t"
                       "vbroadcastsd 208(%%rsi), %%ymm2\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vmovapd %%ymm13, 0(%%rdx)\n\t"
                       "vmovapd %%ymm14, 416(%%rdx)\n\t"
                       "vmovapd %%ymm15, 832(%%rdx)\n\t"
                       "addq $32, %%rdx\n\t"
                       "subq $3712, %%rdi\n\t"
                       "cmpq $52, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $832, %%rdx\n\t"
                       "addq $216, %%rsi\n\t"
                       "subq $416, %%rdi\n\t"
                       "cmpq $9, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 8424;
#endif
}

void libxsmm_m56_n9_k49_ldA56_ldB52_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(const double* A, const double* B, double* C, const double* A_prefetch, const double* B_prefetch, const double* C_prefetch) {
#ifdef __AVX2__
#ifdef __AVX512F__
#pragma message ("LIBXSMM KERNEL COMPILATION WARNING: compiling AVX2 code on AVX512 or newer architecture: " __FILE__)
#endif
  __asm__ __volatile__("movq %0, %%rdi\n\t"
                       "movq %1, %%rsi\n\t"
                       "movq %2, %%rdx\n\t"
                       "movq $0, %%r11\n\t"
                       "33:\n\t"
                       "addq $3, %%r11\n\t"
                       "movq $0, %%r10\n\t"
                       "34:\n\t"
                       "addq $16, %%r10\n\t"
                       "vmovapd 0(%%rdx), %%ymm4\n\t"
                       "vmovapd 32(%%rdx), %%ymm5\n\t"
                       "vmovapd 64(%%rdx), %%ymm6\n\t"
                       "vmovapd 96(%%rdx), %%ymm7\n\t"
                       "vmovapd 448(%%rdx), %%ymm8\n\t"
                       "vmovapd 480(%%rdx), %%ymm9\n\t"
                       "vmovapd 512(%%rdx), %%ymm10\n\t"
                       "vmovapd 544(%%rdx), %%ymm11\n\t"
                       "vmovapd 896(%%rdx), %%ymm12\n\t"
                       "vmovapd 928(%%rdx), %%ymm13\n\t"
                       "vmovapd 960(%%rdx), %%ymm14\n\t"
                       "vmovapd 992(%%rdx), %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 416(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 832(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 416(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 832(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 416(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 832(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 416(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 832(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "cmpq $48, %%r12\n\t"
                       "jl 35b\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 416(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 832(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm4\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm8\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm12\n\t"
                       "vmovapd 32(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm9\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vmovapd 64(%%rdi), %%ymm3\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm6\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm14\n\t"
                       "vmovapd 96(%%rdi), %%ymm3\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm15\n\t"
                       "subq $392, %%rsi\n\t"
                       "vmovapd %%ymm4, 0(%%rdx)\n\t"
                       "vmovapd %%ymm5, 32(%%rdx)\n\t"
                       "vmovapd %%ymm6, 64(%%rdx)\n\t"
                       "vmovapd %%ymm7, 96(%%rdx)\n\t"
                       "vmovapd %%ymm8, 448(%%rdx)\n\t"
                       "vmovapd %%ymm9, 480(%%rdx)\n\t"
                       "vmovapd %%ymm10, 512(%%rdx)\n\t"
                       "vmovapd %%ymm11, 544(%%rdx)\n\t"
                       "vmovapd %%ymm12, 896(%%rdx)\n\t"
                       "vmovapd %%ymm13, 928(%%rdx)\n\t"
                       "vmovapd %%ymm14, 960(%%rdx)\n\t"
                       "vmovapd %%ymm15, 992(%%rdx)\n\t"
                       "addq $128, %%rdx\n\t"
                       "subq $21824, %%rdi\n\t"
                       "cmpq $32, %%r10\n\t"
                       "jl 34b\n\t"
                       "movq $32, %%r10\n\t"
                       "34:\n\t"
                       "addq $12, %%r10\n\t"
                       "vmovapd 0(%%rdx), %%ymm7\n\t"
                       "vmovapd 32(%%rdx), %%ymm8\n\t"
                       "vmovapd 64(%%rdx), %%ymm9\n\t"
                       "vmovapd 448(%%rdx), %%ymm10\n\t"
                       "vmovapd 480(%%rdx), %%ymm11\n\t"
                       "vmovapd 512(%%rdx), %%ymm12\n\t"
                       "vmovapd 896(%%rdx), %%ymm13\n\t"
                       "vmovapd 928(%%rdx), %%ymm14\n\t"
                       "vmovapd 960(%%rdx), %%ymm15\n\t"
                       "movq $0, %%r12\n\t"
                       "35:\n\t"
                       "addq $4, %%r12\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 416(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 832(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 416(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 832(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 416(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 832(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 416(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 832(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "cmpq $48, %%r12\n\t"
                       "jl 35b\n\t"
                       "vbroadcastsd 0(%%rsi), %%ymm0\n\t"
                       "vbroadcastsd 416(%%rsi), %%ymm1\n\t"
                       "vbroadcastsd 832(%%rsi), %%ymm2\n\t"
                       "addq $8, %%rsi\n\t"
                       "vmovapd 0(%%rdi), %%ymm3\n\t"
                       "vmovapd 32(%%rdi), %%ymm4\n\t"
                       "vmovapd 64(%%rdi), %%ymm5\n\t"
                       "vfmadd231pd %%ymm3, %%ymm0, %%ymm7\n\t"
                       "vfmadd231pd %%ymm3, %%ymm1, %%ymm10\n\t"
                       "vfmadd231pd %%ymm3, %%ymm2, %%ymm13\n\t"
                       "vfmadd231pd %%ymm4, %%ymm0, %%ymm8\n\t"
                       "vfmadd231pd %%ymm4, %%ymm1, %%ymm11\n\t"
                       "vfmadd231pd %%ymm4, %%ymm2, %%ymm14\n\t"
                       "addq $448, %%rdi\n\t"
                       "vfmadd231pd %%ymm5, %%ymm0, %%ymm9\n\t"
                       "vfmadd231pd %%ymm5, %%ymm1, %%ymm12\n\t"
                       "vfmadd231pd %%ymm5, %%ymm2, %%ymm15\n\t"
                       "subq $392, %%rsi\n\t"
                       "vmovapd %%ymm7, 0(%%rdx)\n\t"
                       "vmovapd %%ymm8, 32(%%rdx)\n\t"
                       "vmovapd %%ymm9, 64(%%rdx)\n\t"
                       "vmovapd %%ymm10, 448(%%rdx)\n\t"
                       "vmovapd %%ymm11, 480(%%rdx)\n\t"
                       "vmovapd %%ymm12, 512(%%rdx)\n\t"
                       "vmovapd %%ymm13, 896(%%rdx)\n\t"
                       "vmovapd %%ymm14, 928(%%rdx)\n\t"
                       "vmovapd %%ymm15, 960(%%rdx)\n\t"
                       "addq $96, %%rdx\n\t"
                       "subq $21856, %%rdi\n\t"
                       "cmpq $56, %%r10\n\t"
                       "jl 34b\n\t"
                       "addq $896, %%rdx\n\t"
                       "addq $1248, %%rsi\n\t"
                       "subq $448, %%rdi\n\t"
                       "cmpq $9, %%r11\n\t"
                       "jl 33b\n\t"
                       : : "m"(A), "m"(B), "m"(C) : "rdi","rsi","rdx","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15");
#else
#pragma message ("LIBXSMM KERNEL COMPILATION ERROR in: " __FILE__)
#error No kernel was compiled, lacking support for current architecture?
#endif

#ifndef NDEBUG
#ifdef _OPENMP
#pragma omp atomic
#endif
libxsmm_num_total_flops += 49392;
#endif
}

