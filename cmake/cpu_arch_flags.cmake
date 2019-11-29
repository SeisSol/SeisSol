function(get_arch_flags architecture compiler)
    # Westmere cpu architecture
    if ("${ARCH}" STREQUAL "wsm")
        set(CPU_ARCH_FLAGS "-msse3" PARENT_SCOPE)
    
    # Sandy Bridge cpu architecture
    elseif ("${ARCH}" STREQUAL "snb")
        set(CPU_ARCH_FLAGS "-mavx" PARENT_SCOPE)
    
    # Haswell cpu architecture
    elseif ("${ARCH}" STREQUAL "hsw")
        if (compiler STREQUAL "Intel")
            set(CPU_ARCH_FLAGS "-xCORE-AVX2" "-fma" PARENT_SCOPE)
        elseif(compiler STREQUAL "GNU")
            set(CPU_ARCH_FLAGS "-mavx2" "-mfma" PARENT_SCOPE)
        endif()

    # Knights Corner (Xeon Phi)
    elseif ("${ARCH}" STREQUAL "knc")
        set(CPU_ARCH_FLAGS "-mmic" "-fma" PARENT_SCOPE)
    
    # Knight Landing (Xeon Phi)
    elseif ("${ARCH}" STREQUAL "knl")
        if (compiler STREQUAL "Intel")
            set(CPU_ARCH_FLAGS "-xMIC-AVX512" "-fma" PARENT_SCOPE)
        elseif(compiler STREQUAL "GNU")
            set(CPU_ARCH_FLAGS "-mavx512f" "-mavx512cd" "-mavx512pf" "-mavx512er" "-mfma" PARENT_SCOPE)
        endif()
    
    # Skylake cpu architecture
    elseif ("${ARCH}" STREQUAL "skx")
        if (compiler STREQUAL "Intel")
            set(CPU_ARCH_FLAGS "-xMIC-AVX512" "-fma" PARENT_SCOPE)
        elseif(compiler STREQUAL "GNU")
            set(CPU_ARCH_FLAGS "-march=skylake-avx512" PARENT_SCOPE)
        endif()
    endif()

endfunction()