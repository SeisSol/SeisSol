/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 **/

#include "LTS.h"
                                                 // Ghost,  Copy, Interior
bool const LTS::Available[][LTS::NUM_LAYERS] =  { { false,  true,   true  },  // Dofs
                                                  { true,   true,   true  },  // Buffers
                                                  { true,   true,   true  },  // Derivatives
                                                  { false,  true,   true  },  // FaceNeighbors
                                                  { false,  true,   true  },  // LocalIntegration
                                                  { false,  true,   true  },  // NeighboringIntegration
                                                  { false,  true,   true  },  // Material
#ifdef USE_PLASTICITY
                                                  { false,  true,   true  },  // Plasticity
                                                  { false,  true,   true  },  // Energy
                                                  { false,  true,   true  }   // PStrain
#else
                                                  { false,  false,  false },  // Plasticity
                                                  { false,  false,  false },  // Energy
                                                  { false,  false,  false }   // PStrain
#endif
                                                };
size_t const LTS::VarAlignment[] =  { PAGESIZE_HEAP,  // Dofs
                                      1,              // Buffers
                                      1,              // Derivatives
                                      1,              // FaceNeighbors
                                      1,              // LocalIntegration
                                      1,              // NeighboringIntegration
                                      1,              // Material
                                      1,              // Plasticity
                                      PAGESIZE_HEAP,  // Energy
                                      PAGESIZE_HEAP   // PStrain
                                    };
/// \todo Memkind
enum seissol::memory::Memkind const LTS::VarMemkind[] = { seissol::memory::Standard,  // Dofs
                                                          seissol::memory::Standard,  // Buffers
                                                          seissol::memory::Standard,  // Derivatives
                                                          seissol::memory::Standard,  // FaceNeighbors
                                                          seissol::memory::Standard,  // LocalIntegration
                                                          seissol::memory::Standard,  // NeighboringIntegration
                                                          seissol::memory::Standard,  // Material
                                                          seissol::memory::Standard,  // Plasticity
                                                          seissol::memory::Standard,  // Energy
                                                          seissol::memory::Standard   // PStrain
                                                        };
size_t const LTS::BucketAlignment[] = { PAGESIZE_HEAP }; // buffersderivatives
/// \todo Memkind
enum seissol::memory::Memkind const LTS::BucketMemkind[] = { seissol::memory::Standard }; // buffersderivatives
