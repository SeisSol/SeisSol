/*
Copyright (c) 2015, Intel Corporation

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

typedef struct seissol_flops {
  long long d_nonZeroFlops;
  long long d_hardwareFlops;
} seissol_flops;

seissol_flops flops_localWithoutAder_actual(unsigned int i_timesteps) {
  seissol_flops ret;
  ret.d_nonZeroFlops = 0.0;
  ret.d_hardwareFlops = 0.0;

  auto&                 layer           = m_ltsTree.child(0).child<Interior>();
  unsigned              nrOfCells       = layer.getNumberOfCells();
  CellLocalInformation* cellInformation = layer.var(m_lts.cellInformation);
  for (unsigned cell = 0; cell < nrOfCells; ++cell) {
    unsigned int l_nonZeroFlops, l_hardwareFlops;
    m_localKernel.flopsIntegral(cellInformation[cell].faceTypes, l_nonZeroFlops, l_hardwareFlops);    
    ret.d_nonZeroFlops  += l_nonZeroFlops;
    ret.d_hardwareFlops += l_hardwareFlops;
  }

  ret.d_nonZeroFlops *= i_timesteps;
  ret.d_hardwareFlops *= i_timesteps;
  
  return ret;
}

seissol_flops flops_ader_actual(unsigned int i_timesteps) {
  seissol_flops ret;
  ret.d_nonZeroFlops = 0.0;
  ret.d_hardwareFlops = 0.0;
  
  // iterate over cells
  unsigned nrOfCells = m_ltsTree.child(0).child<Interior>().getNumberOfCells();
  for( unsigned int l_cell = 0; l_cell < nrOfCells; l_cell++ ) {
    unsigned int l_nonZeroFlops, l_hardwareFlops;
    // get flops
    m_timeKernel.flopsAder( l_nonZeroFlops, l_hardwareFlops );
    ret.d_nonZeroFlops  += l_nonZeroFlops;
    ret.d_hardwareFlops += l_hardwareFlops;
  }

  ret.d_nonZeroFlops *= i_timesteps;
  ret.d_hardwareFlops *= i_timesteps;

  return ret;
}

seissol_flops flops_neigh_actual(unsigned int i_timesteps) {
  seissol_flops ret;
  ret.d_nonZeroFlops = 0.0;
  ret.d_hardwareFlops = 0.0;
  
  // iterate over cells
  auto&                 layer           = m_ltsTree.child(0).child<Interior>();
  unsigned              nrOfCells       = layer.getNumberOfCells();
  CellLocalInformation* cellInformation = layer.var(m_lts.cellInformation);
  CellDRMapping        (*drMapping)[4]  = layer.var(m_lts.drMapping);
  for( unsigned int l_cell = 0; l_cell < nrOfCells; l_cell++ ) {
    unsigned int l_nonZeroFlops, l_hardwareFlops;
    long long l_drNonZeroFlops, l_drHardwareFlops;
    // get flops
    m_neighborKernel.flopsNeighborsIntegral( cellInformation[l_cell].faceTypes, cellInformation[l_cell].faceRelations, drMapping[l_cell], l_nonZeroFlops, l_hardwareFlops, l_drNonZeroFlops, l_drHardwareFlops );
    ret.d_nonZeroFlops  += l_nonZeroFlops + l_drNonZeroFlops;
    ret.d_hardwareFlops += l_hardwareFlops + l_drHardwareFlops;
  }

  ret.d_nonZeroFlops *= i_timesteps;
  ret.d_hardwareFlops *= i_timesteps;

  return ret;
}

seissol_flops flops_drgod_actual(unsigned int i_timesteps) {
  seissol_flops ret;
  ret.d_nonZeroFlops = 0.0;
  ret.d_hardwareFlops = 0.0;
  
  // iterate over cells
  seissol::initializers::Layer& interior = m_dynRupTree.child(0).child<Interior>();
  DRFaceInformation* faceInformation = interior.var(m_dynRup.faceInformation);
  for (unsigned face = 0; face < interior.getNumberOfCells(); ++face) {
    long long l_drNonZeroFlops, l_drHardwareFlops;
    m_dynRupKernel.flopsGodunovState(faceInformation[face], l_drNonZeroFlops, l_drHardwareFlops);
    ret.d_nonZeroFlops  += l_drNonZeroFlops;
    ret.d_hardwareFlops += l_drHardwareFlops;
  }

  ret.d_nonZeroFlops *= i_timesteps;
  ret.d_hardwareFlops *= i_timesteps;

  return ret;
}

seissol_flops flops_local_actual(unsigned int i_timesteps) {
  seissol_flops ret;
  seissol_flops tmp;
  
  tmp = flops_ader_actual(i_timesteps);
  ret.d_nonZeroFlops = tmp.d_nonZeroFlops;
  ret.d_hardwareFlops = tmp.d_hardwareFlops;

  tmp = flops_localWithoutAder_actual(i_timesteps);
  ret.d_nonZeroFlops += tmp.d_nonZeroFlops;
  ret.d_hardwareFlops += tmp.d_hardwareFlops;

  return ret;
}

seissol_flops flops_all_actual(unsigned int i_timesteps) {
  seissol_flops ret;
  seissol_flops tmp;
  
  tmp = flops_local_actual(i_timesteps);
  ret.d_nonZeroFlops = tmp.d_nonZeroFlops;
  ret.d_hardwareFlops = tmp.d_hardwareFlops;

  tmp = flops_neigh_actual(i_timesteps);
  ret.d_nonZeroFlops += tmp.d_nonZeroFlops;
  ret.d_hardwareFlops += tmp.d_hardwareFlops;

  return ret;
}


