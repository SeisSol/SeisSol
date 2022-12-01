function(get_arch_flags architecture compiler)
    set(HAS_REDZONE ON PARENT_SCOPE)

    # Westmere cpu architecture
    if ("${HOST_ARCH}" STREQUAL "wsm")
        set(CPU_ARCH_FLAGS "-msse3" PARENT_SCOPE)
    
    # Sandy Bridge cpu architecture
    elseif ("${HOST_ARCH}" STREQUAL "snb")
        set(CPU_ARCH_FLAGS "-mavx" PARENT_SCOPE)
    
    # Haswell cpu architecture
    elseif ("${HOST_ARCH}" STREQUAL "hsw")
        if (compiler STREQUAL "Intel")
            set(CPU_ARCH_FLAGS "-xCORE-AVX2" "-fma" PARENT_SCOPE)
        elseif(compiler MATCHES "GNU|Clang|IntelLLVM")
            set(CPU_ARCH_FLAGS "-mavx2" "-mfma" PARENT_SCOPE)
        endif()

    # Knights Corner (Xeon Phi)
    elseif ("${HOST_ARCH}" STREQUAL "knc")
        set(CPU_ARCH_FLAGS "-mmic" "-fma" PARENT_SCOPE)
    
    # Knight Landing (Xeon Phi)
    elseif ("${HOST_ARCH}" STREQUAL "knl")
        if (compiler STREQUAL "Intel")
            set(CPU_ARCH_FLAGS "-xMIC-AVX512" "-fma" PARENT_SCOPE)
        elseif(compiler MATCHES "GNU|Clang|IntelLLVM")
            set(CPU_ARCH_FLAGS "-mavx512f" "-mavx512cd" "-mavx512pf" "-mavx512er" "-mfma" PARENT_SCOPE)
        endif()
    
    # Skylake cpu architecture
    elseif ("${HOST_ARCH}" STREQUAL "skx")
        if (compiler STREQUAL "Intel")
            set(CPU_ARCH_FLAGS "-xCORE-AVX512" "-fma" PARENT_SCOPE)
        elseif(compiler MATCHES "GNU|Clang|IntelLLVM")
            set(CPU_ARCH_FLAGS "-march=skylake-avx512" PARENT_SCOPE)
        endif()

    elseif ("${HOST_ARCH}" STREQUAL "thunderx2t99")
        set(HAS_REDZONE OFF PARENT_SCOPE)
        if (compiler STREQUAL "Intel")

        elseif(compiler MATCHES "GNU|Clang")
	    # Note: mcpu/march/mtune are weird on arm, see:
	    # https://community.arm.com/developer/tools-software/tools/b/tools-software-ides-blog/posts/compiler-flags-across-architectures-march-mtune-and-mcpu
            set(CPU_ARCH_FLAGS "-mcpu=thunderx2t99" PARENT_SCOPE)
        endif()
    # AMD Rome/ Epyc 2nd Gen
    elseif ("${HOST_ARCH}" STREQUAL "rome")
        if (compiler STREQUAL "Intel")
            set(CPU_ARCH_FLAGS "-march=core-avx2" "-fma" PARENT_SCOPE)
        elseif(compiler MATCHES "GNU|Clang|IntelLLVM")
            set(CPU_ARCH_FLAGS "-march=znver2" "-mtune=znver2" PARENT_SCOPE)
        endif()

    # IBM power 9
    elseif ("${HOST_ARCH}" STREQUAL "power9")
        if (compiler MATCHES "GNU|Clang")
            set(HAS_REDZONE OFF PARENT_SCOPE)
            set(CPU_ARCH_FLAGS "-mtune=power9" PARENT_SCOPE)
        endif()
    endif()

    if (compiler MATCHES "NVHPC|PGI")
        #NOTE: PGI-based compiler does not have `-mno-red-zone` flag
        set(HAS_REDZONE OFF PARENT_SCOPE)
    endif()

endfunction()
