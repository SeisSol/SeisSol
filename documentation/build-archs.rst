..
  SPDX-FileCopyrightText: 2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _build_archs:

Build architectures
===================

Currently, SeisSol needs information about the host architecture on which the code is going to run.

Besides setting the necessary compiler tuning variables (usually corresponding to ``-march=TARGET_ARCH -mtune=TARGET_ARCH``),
it also sets the code generators.

CPU architectures
~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 20 40 20 20 20
   :header-rows: 1

   * - ``HOST_ARCH``
     - Architecture
     - Notes
     - CPUs (examples)
     - Usage
   * - ``noarch``
     - No architecture-specific optimizations
     - Generates plain x86-64 instructions, without SIMD instructions like SSE/AVX/AMX etc.
     -
     -
   * - ``wsm``
     - Intel Nehalem/Westmere architecture
     - Generates SSE instructions (up to SSE 3).
     - Intel Xeon v1/v0, Intel Core i3/i5/i7 ??? and 1???
     -
   * - ``snb``
     - Intel Sandy Bridge architecture
     - Generates AVX instructions.
     - Intel Xeon v2, Intel Core i3/i5/i7 2??? and 3???
     -
   * - ``hsw``
     - Intel Haswell
     - Generates AVX2 instructions.
     - Intel Xeon v3, Intel Core i3/i5/i7/i9 4??? to 14??? (mostly), AMD Zen 1 to 3
     - Older Intel CPU clusters, clusters with AMD CPUs (up to 2023)
   * - ``skx``
     - Intel Skylake-X (including Skylake-SP)
     - Generates AVX-512{F,CD,BW,DQ,VL} instructions. (NOTE: Skylake desktop processors are NOT included here, unless they contain an "X" in their name, such as e.g. i9 7800X)
     - Intel Xeon v4 and onward (i.e. including the "metal"-branded Xeons), some Intel Core i9 models (check the Intel database), AMD Zen 4
     - Most CPU clusters, e.g. SuperMUC NG (Phase 1), Frontera
   * - ``knc``
     - Intel Knight's Corner (Xeon Phi coprocessor)
     - Generates Knight's Corner-specific instructions.
     - Intel Xeon Phi coprocessor
     - (not known anymore)
   * - ``knl``
     - Intel Knight's Landing (Xeon Phi, optionally as coprocessor)
     - Generates AVX-512{F,CD,PF,ER} instructions.
     - Intel Xeon Phi coprocessor, as well as
     - LRZ CoolMUC 3
   * - ``naples``
     - AMD Zen 1
     - Generates AVX2 instructions. For the libxsmm kernel generator, it is deemed equivalent to ``hsw``.
     - Ryzen 1xxx series
     -
   * - ``rome``
     - AMD Zen 2
     - Generates AVX2 instructions. For the libxsmm kernel generator, it is deemed equivalent to ``hsw``.
     - Ryzen 3??? series, 7?2?, 8?2? series
     - LUMI (CPU partition)
   * - ``milan``
     - AMD Zen 3
     - Generates AVX2 instructions. For the libxsmm kernel generator, it is deemed equivalent to ``hsw``.
     - Ryzen 5??? series, 7?3?, 8?3? series
     - LUMI (GPU partition), Frontier (GPU partition)
   * - ``bergamo``
     - AMD Zen 4
     - Generates AVX512 instructions. For the libxsmm kernel generator, it is deemed equivalent to ``skx``.
     - Ryzen 7?4? series, MI300A
     -
   * - ``power9``
     - IBM PowerPC 9
     -
     -
     -
   * - ``thunderx2t99``
     - ARM ThunderX2 (ARM NEON)
     - ARM ThunderX2
     - Isambard 2
     -
   * - ``a64fx``
     - Fujitsu A64FX (ARM SVE, 512 bits)
     - Fujitsu A64FX
     - Fugaku
     -
   * - ``neon``
     - Dummy target for AARCH64 (with NEON)
     -
     -
     -
   * - ``sve128``
     - Dummy target for AARCH64, ARM SVE with 128 bits length
     -
     - Needed e.g. for the Neoverse V2 CPU
     -
   * - ``sve256``
     - Dummy target for AARCH64, ARM SVE with 256 bits length
     -
     -
     -
   * - ``sve512``
     - Dummy target for AARCH64, ARM SVE with 512 bits length
     -
     -
     -
   * - ``sve1024``
     - Dummy target for AARCH64, ARM SVE with 1024 bits length
     -
     -
     -
   * - ``sve2048``
     - Dummy target for AARCH64, ARM SVE with 2048 bits length
     -
     -
     -
   * - ``apple-m1``
     - Apple M1 CPU
     -
     -
     -
   * - ``apple-m2``
     - Apple M2 CPU
     -
     -
     -

GPU architectures
~~~~~~~~~~~~~~~~~

For GPUs, SeisSol supports two types of memory management on GPUs.

* **split**: separate host and device buffers which are synchronized for e.g. IO
* **unified**: combined host-device buffers which are transferred by the system as needed

As default, SeisSol will use unified host-device buffers by default on all systems where the CPU can freely access
GPU memory, i.e. the Nvidia Superchips (e.g. GH100, e.g. GH200) and the AMD APUs (e.g. MI300A).

In all other cases, split host-device buffers will be used as default.

The following architectures are supported:

.. list-table::
   :widths: 20 40 20 20 20
   :header-rows: 1

   * - ``DEVICE_ARCH``
     - ``DEVICE_BACKEND``
     - Architecture
     - GPUs (examples)
     - Memory default
   * - ``sm_60``
     - ``cuda``
     - Nvidia Pascal
     - Nvidia P100
     - split
   * - ``sm_61``
     - ``cuda``
     - Nvidia Pascal
     - Nvidia Geforce 1000 series, Quadro P series
     - split
   * - ``sm_70``
     - ``cuda``
     - Nvidia Volta
     - Nvidia V100
     - split
   * - ``sm_75``
     - ``cuda``
     - Nvidia Turing
     - Nvidia Geforce 2000 series, Quadro RTX series
     - split
   * - ``sm_80``
     - ``cuda``
     - Nvidia Ampere
     - Nvidia A100
     - split
   * - ``sm_86``
     - ``cuda``
     - Nvidia Ampere
     - Nvidia Geforce 3000 series, Quadro RTX A series
     - split
   * - ``sm_89``
     - ``cuda``
     - Nvidia Lovelace
     - Nvidia Geforce 4000 series, Quadro RTX Ada series
     - split
   * - ``sm_90``
     - ``cuda``
     - Nvidia Hopper
     - Nvidia H100, H200
     - split; unified on GH superchip
   * - ``sm_100``
     - ``cuda``
     - Nvidia Blackwell
     - Nvidia B100, B200
     - split; unified on GB superchip
   * - ``gfx900``
     - ``hip``
     - AMD GCN 5 (Vega)
     - AMD Instinct MI25, Radeon RX Vega 56, Radeon RX Vega 64
     - split [#xnack1]_
   * - ``gfx906``
     - ``hip``
     - AMD GCN 5 (Vega)
     - AMD Instinct MI50, Radeon VII
     - split [#xnack1]_
   * - ``gfx908``
     - ``hip``
     - AMD CDNA 1
     - AMD Instinct MI100X
     - split [#xnack1]_
   * - ``gfx90a``
     - ``hip``
     - AMD CDNA 2
     - AMD Instinct MI210, MI250X
     - split [#xnack1]_
   * - ``gfx942``
     - ``hip``
     - AMD CDNA 3
     - AMD Instinct MI300A, MI300X
     - split [#xnack1]_; unified on MI300A
   * - ``gfx1010``
     - ``hip``
     - AMD RDNA 1
     - AMD Radeon 5000 series
     - split [#xnack1]_
   * - ``gfx1030``
     - ``hip``
     - AMD RDNA 2
     - AMD Radeon 6000 series
     - split [#xnack2]_
   * - ``gfx1100``
     - ``hip``
     - AMD RDNA 3
     - AMD Radeon 7000 series
     - split [#xnack2]_
   * - ``pvc``
     - ``oneapi``
     - Intel Ponte Vecchio
     - Intel Data Center Max 1550
     - split

Sources:

* https://en.wikipedia.org/wiki/CUDA#GPUs_supported
* https://llvm.org/docs/AMDGPUUsage.html
* https://intel.github.io/llvm-docs/UsersManual.html

About AMD GPUs: for unified memory to perform well, you will need to set ``HSA_XNACK=1``.

For unsupported AMD GPU architectures (e.g. ``gfx90c``), you can proceed as follows:

* compile for a compatible GPU architecture. In the case of ``gfx90c``, your best choice is ``gfx900`` (or ``gfx906``).
* run SeisSol with specifying the environment variable ``HSA_OVERRIDE_GFX_VERSION`` in accordance to the architecture you compiled against in the previous step. That is, you need to convert ``gfxAABC`` to a version of the form ``AA.B.C``. E.g., if you compiled for ``gfx906``, you will need to set ``HSA_OVERRIDE_GFX_VERSION=9.0.6``. Letters become numbers, akin to the hexadecimal notation, i.e. ``gfx90a`` becomes 9.0.10.

.. [#xnack1] For managed memory support to perform well, you will need to set ``HSA_XNACK=1`` as environment variable.
.. [#xnack2] Do not support managed memory.
