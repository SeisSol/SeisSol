#include <cstdio>

#include <Kernels/Time.h>
#include <Kernels/Local.h>
#include <Kernels/Neighbor.h>
#include <Kernels/DynamicRupture.h>

int main()
{
  /// ADER-DG Flops
  long long nonZeroFlops = 0, hardwareFlops = 0;
  long long neighborNonZeroFlops = 0, neighborHardwareFlops = 0;
  unsigned kernelNonZeroFlops, kernelHardwareFlops;
  long long kernelDRNonZeroFlops, kernelDRHardwareFlops;

  seissol::kernels::Time      timeKernel;
  seissol::kernels::Local     localKernel;
  seissol::kernels::Neighbor  neighborKernel;
  
  faceType faceTypes[] = {regular, regular, regular, regular};
  
  timeKernel.flopsAder( kernelNonZeroFlops, kernelHardwareFlops );
  nonZeroFlops += kernelNonZeroFlops;
  hardwareFlops += kernelHardwareFlops;

  localKernel.flopsIntegral(faceTypes, kernelNonZeroFlops, kernelHardwareFlops);
  nonZeroFlops += kernelNonZeroFlops;
  hardwareFlops += kernelHardwareFlops;
  
  int faceRelations[4][2];
  
  for (int neighborSide = 0; neighborSide < 4; ++neighborSide) {
    for (int sideOrientation = 0; sideOrientation < 3; ++sideOrientation) {
      for (int side = 0; side < 4; ++side) {
        faceRelations[side][0] = neighborSide;
        faceRelations[side][1] = sideOrientation;
      }
      neighborKernel.flopsNeighborsIntegral( faceTypes, faceRelations, kernelNonZeroFlops, kernelHardwareFlops, kernelDRNonZeroFlops, kernelDRHardwareFlops );
      neighborNonZeroFlops += kernelNonZeroFlops;
      neighborHardwareFlops += kernelHardwareFlops;
    }
  }
  
  printf("Element non-zero flops (average): %.1lf\n", nonZeroFlops + neighborNonZeroFlops / 12.0);
  printf("Element hardware flops (average): %.1lf\n", hardwareFlops + neighborHardwareFlops / 12.0);

  
  /// Dynamic rupture flops

  long long drNonZeroFlops = 0, drHardwareFlops = 0;
  long long neighborDRNonZeroFlops = 0, neighborDRHardwareFlops = 0;
  
  seissol::kernels::DynamicRupture dynRupKernel;

  faceType faceTypesDR[] = {dynamicRupture, dynamicRupture, dynamicRupture, dynamicRupture};
  
  for (int neighborSide = 0; neighborSide < 4; ++neighborSide) {
    for (int sideOrientation = 0; sideOrientation < 3; ++sideOrientation) {
      for (int side = 0; side < 4; ++side) {
        DRFaceInformation faceInfo;
        faceInfo.plusSide = side;
        faceInfo.minusSide = neighborSide;
        faceInfo.faceRelation = sideOrientation;
        
        dynRupKernel.flopsGodunovState(faceInfo, kernelDRNonZeroFlops, kernelDRHardwareFlops);
        drNonZeroFlops += kernelDRNonZeroFlops;
        drHardwareFlops += kernelDRHardwareFlops;
      }
    }
  }
  
  for (int neighborSide = 0; neighborSide < 4; ++neighborSide) {
    for (int sideOrientation = 0; sideOrientation < 3; ++sideOrientation) {
      for (int face = 0; face < 4; ++face) {
        faceRelations[face][0] = neighborSide;
        faceRelations[face][1] = sideOrientation;
      }
      neighborKernel.flopsNeighborsIntegral( faceTypesDR, faceRelations, kernelNonZeroFlops, kernelHardwareFlops, kernelDRNonZeroFlops, kernelDRHardwareFlops );
      neighborDRNonZeroFlops += kernelDRNonZeroFlops;
      neighborDRHardwareFlops += kernelDRHardwareFlops;
    }
  }  
  
  printf("Dynamic rupture face non-zero flops (average, w/o friction law): %.1lf\n", drNonZeroFlops / 48.0 + neighborDRNonZeroFlops / 12.0 / 4.0);
  printf("Dynamic rupture face hardware flops (average, w/o friction law): %.1lf\n", drHardwareFlops / 48.0 + neighborDRHardwareFlops / 12.0 / 4.0);
  
  return 0;
}
