#include <cstdio>

#include <Kernels/Time.h>
#include <Kernels/Local.h>
#include <Kernels/Neighbor.h>
#include <Kernels/DynamicRupture.h>
#include <Kernels/Plasticity.h>

int main()
{
  /// ADER-DG Flops
  long long nonZeroFlops = 0, hardwareFlops = 0;
  long long avgNeighborNonZeroFlops = 0, avgNeighborHardwareFlops = 0;
  long long minNeighborNonZeroFlops = std::numeric_limits<long long>::max(), minNeighborHardwareFlops = std::numeric_limits<long long>::max();
  long long maxNeighborNonZeroFlops = 0, maxNeighborHardwareFlops = 0;
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
  CellDRMapping dummyDRMapping[4];

  unsigned numConfs = 12*12*12*12;
  for (unsigned conf = 0; conf < 12*12*12*12; ++conf) {
    unsigned m = 1;
    for (int side = 0; side < 4; ++side) {
      unsigned confSide = (conf % (12*m)) / m;
      faceRelations[side][0] = confSide % 4;
      faceRelations[side][1] = confSide / 4;
      m *= 12;
    }
    neighborKernel.flopsNeighborsIntegral( faceTypes, faceRelations, dummyDRMapping, kernelNonZeroFlops, kernelHardwareFlops, kernelDRNonZeroFlops, kernelDRHardwareFlops );
    minNeighborNonZeroFlops = std::min(minNeighborNonZeroFlops, static_cast<long long>(kernelNonZeroFlops));
    maxNeighborNonZeroFlops = std::max(maxNeighborNonZeroFlops, static_cast<long long>(kernelNonZeroFlops));
    minNeighborHardwareFlops = std::min(minNeighborHardwareFlops, static_cast<long long>(kernelHardwareFlops));
    maxNeighborHardwareFlops = std::max(maxNeighborHardwareFlops, static_cast<long long>(kernelHardwareFlops));
    avgNeighborNonZeroFlops += kernelNonZeroFlops;
    avgNeighborHardwareFlops += kernelHardwareFlops;
  }
  
  double avgNonZeroFlops = nonZeroFlops + avgNeighborNonZeroFlops / static_cast<double>(numConfs);
  double avgHardwareFlops = hardwareFlops + avgNeighborHardwareFlops / static_cast<double>(numConfs);
  printf("Element non-zero flops (average): %.1lf\n", avgNonZeroFlops);
  printf("Element hardware flops (average): %.1lf\n", avgHardwareFlops);
  printf("Element non-zero flops (min): %lld\n", nonZeroFlops + minNeighborNonZeroFlops);
  printf("Element hardware flops (min): %lld\n", hardwareFlops + minNeighborHardwareFlops);
  printf("Element non-zero flops (max): %lld\n", nonZeroFlops + maxNeighborNonZeroFlops);
  printf("Element hardware flops (max): %lld\n", hardwareFlops + maxNeighborHardwareFlops);
  printf("Max non-zero deviation from average: %.2lf %\n", 100.0 * std::max((nonZeroFlops + maxNeighborNonZeroFlops) / avgNonZeroFlops - 1.0, 1.0 - (nonZeroFlops + minNeighborNonZeroFlops) / avgNonZeroFlops));
  printf("Max hardware deviation from average: %.2lf %\n", 100.0 * std::max((hardwareFlops + maxNeighborHardwareFlops) / avgHardwareFlops - 1.0, 1.0 - (hardwareFlops + minNeighborHardwareFlops) / avgHardwareFlops));


  
  /// Dynamic rupture flops

  long long drNonZeroFlops = 0, drHardwareFlops = 0;
  long long neighborDRNonZeroFlops = 0, neighborDRHardwareFlops = 0;
  
  seissol::kernels::DynamicRupture dynRupKernel;

  faceType faceTypesDR[] = {dynamicRupture, dynamicRupture, dynamicRupture, dynamicRupture};
  CellDRMapping drMapping[4];
  
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
    for (int faceRelation = 0; faceRelation < 4; ++faceRelation) {
      for (int face = 0; face < 4; ++face) {
        drMapping[face].side = neighborSide;
        drMapping[face].faceRelation = faceRelation;
      }
      neighborKernel.flopsNeighborsIntegral( faceTypesDR, faceRelations, drMapping, kernelNonZeroFlops, kernelHardwareFlops, kernelDRNonZeroFlops, kernelDRHardwareFlops );
      neighborDRNonZeroFlops += kernelDRNonZeroFlops;
      neighborDRHardwareFlops += kernelDRHardwareFlops;
    }
  }  
  
  printf("Dynamic rupture face non-zero flops (average, w/o friction law): %.1lf\n", drNonZeroFlops / 48.0 + neighborDRNonZeroFlops / 16.0 / 4.0);
  printf("Dynamic rupture face hardware flops (average, w/o friction law): %.1lf\n", drHardwareFlops / 48.0 + neighborDRHardwareFlops / 16.0 / 4.0);
  

  /// Plasticity flops

  long long PlNonZeroFlopsCheck = 0, PlHardwareFlopsCheck = 0;
  long long PlNonZeroFlopsYield = 0, PlHardwareFlopsYield = 0;
  long long kernelNonZeroFlopsCheck, kernelHardwareFlopsCheck;
  long long kernelNonZeroFlopsYield, kernelHardwareFlopsYield;

  seissol::kernels::Plasticity::flopsPlasticity(kernelNonZeroFlopsCheck, kernelHardwareFlopsCheck, kernelNonZeroFlopsYield, kernelHardwareFlopsYield);

  PlNonZeroFlopsCheck += kernelNonZeroFlopsCheck;
  PlHardwareFlopsCheck += kernelHardwareFlopsCheck;
  PlNonZeroFlopsYield += kernelNonZeroFlopsYield;
  PlHardwareFlopsYield += kernelHardwareFlopsYield;

  printf("Plasticity non-zero min flops : %lld\n", PlNonZeroFlopsCheck );
  printf("Plasticity hardware min flops : %lld\n", PlHardwareFlopsCheck );
  printf("Plasticity non-zero max flops : %lld\n", PlNonZeroFlopsCheck + PlNonZeroFlopsYield );
  printf("Plasticity hardware max flops : %lld\n", PlHardwareFlopsCheck+ PlHardwareFlopsYield );

  printf("Plasticity non-zero min vs elastic average: %.2lf %\n", 100.0 * PlNonZeroFlopsCheck / avgNonZeroFlops);
  printf("Plasticity hardware min vs elastic average: %.2lf %\n", 100.0 * PlHardwareFlopsCheck / avgHardwareFlops);
  printf("Plasticity non-zero max vs elastic average: %.2lf %\n", 100.0 * (PlNonZeroFlopsCheck + PlNonZeroFlopsYield) / avgNonZeroFlops);
  printf("Plasticity hardware max vs elastic average: %.2lf %\n", 100.0 * (PlHardwareFlopsCheck + PlHardwareFlopsYield) / avgHardwareFlops);

  return 0;
}
