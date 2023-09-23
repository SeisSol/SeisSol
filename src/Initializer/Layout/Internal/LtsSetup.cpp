#include "LtsSetup.hpp"
#include <Initializer/Layout/Memory.hpp>
#include <Initializer/typedefs.hpp>
#include <array>

namespace {


/**
 * Gets the lts setup in relation to the four face neighbors.
 *   Remark: Remember to perform the required normalization step.
 *
 * -------------------------------------------------------------------------------
 *
 *  0 in one of the first four bits: Face neighboring data are buffers.
 *  1 in one of the first four bits: Face neighboring data are derivatives.
 *
 *     Example 1:
 *     [           12 rem. bits               | buf/der bits ]
 *     [  -  -  -  -  -  -  -  -  -  -  -  -  |  0  1  1  0  ]
 *     [ 15 14 13 12 11 10  9  8  7  6  5  4  |  3  2  1  0  ]
 *  In Example 1 the data for face neighbors 0 and 3 are buffers and for 1 and 2 derivatives.
 *
 *  0 in one of bits 4 - 7: No global time stepping
 *  1 in one of bits 4 - 7: The  current cell has a global time stepping relation with the face neighbor.
 *
 *     Example 2:
 *     [       8 rem. bits       |   GTS bits  | buf/der bits ]
 *     [  -  -  -  -  -  -  -  - | 0  0  1  1  |  0  1  1  0  ]
 *     [ 15 14 13 12 11 10  9  8 | 7  6  5  4  |  3  2  1  0  ]
 *  In Example 2 the data of face neighbors 0 and 3 are buffers, 1 and 2 deliver derivatives
 *  Face neighbor 0 has a GTS-relation and this cell works directly on the delivered buffer.
 *  Face neighbor 1 has a GTS-relation, but delivers derivatives -> The derivatives have to translated to time integrated DOFs first.
 *  Face neighbor 2 has a LTS-relation and receives derivatives from its neighbor -> The derivates have to be used for a partial time integration.
 *  Face neighbor 3 has a LTS-relation and can operate on the buffers directly.
 *
 * -------------------------------------------------------------------------------
 *
 *  1 in the eigth bit: the cell is required to work on time integration buffers.
 *  1 in the nineth bit: the cell is required to compute time derivatives.
 *
 *     Example 3:
 *     [     remaining     | der. buf. |       first 8 bits       ]
 *     [  -  -  -  -  -  - |  0    1   |  -  -  -  -  -  -  -  -  ]
 *     [ 15 14 13 12 11 10 |  9    8   |  7  6  5  4  3  2  1  0  ]
 *  In Example 3 only a buffer is stored as for example in global time stepping.
 *
 *     Example 4:
 *     [     remaining     | der. buf. |       first 8 bits       ]
 *     [  -  -  -  -  -  - |  1    1   |  -  -  -  -  -  -  -  -  ]
 *     [ 15 14 13 12 11 10 |  9    8   |  7  6  5  4  3  2  1  0  ]
 *  In Example 4 both (buffer+derivative) is stored.
 *
 * -------------------------------------------------------------------------------
 *
 *  1 in the tenth bit: the cell local buffer is a LTS buffer (reset on request only).
 *
 *     Example 5:
 *     [   remaining    | LTS buf. |          first 10 bits         ]
 *     [  -  -  -  -  - |     1    |  -  1  -  -  -  -  -  -  -  -  ]
 *     [ 15 14 13 12 11 |    10    |  9  8  7  6  5  4  3  2  1  0  ]
 *  In Example 5 the buffer is a LTS buffer (reset on request only). GTS buffers are updated in every time step.
 *
 *  Proposed extension (WIP): bits 11-14 specify if we need to project to a nodal basis on the faces
 **/
    unsigned short determineLtsSetup(const CellLocalInformation& ownPrimary, const SecondaryCellLocalInformation& ownSecondary, const std::array<std::reference_wrapper<CellLocalInformation>, 4>& primary, const std::array<std::reference_wrapper<SecondaryCellLocalInformation>, 4>& secondary, bool copyLayer) {
// reset the LTS setup
  unsigned short ltsSetup = 0;

  // iterate over the faces
  for( unsigned int face = 0; face < 4; face++ ) {
    // continue for boundary conditions
    if (ownPrimary.faceTypes[face] == FaceType::outflow) {
      continue;
    }
    // fake neighbors are GTS
    else if(ownPrimary.faceTypes[face] == FaceType::freeSurface ||
	    ownPrimary.faceTypes[face] == FaceType::freeSurfaceGravity ||
	    ownPrimary.faceTypes[face] == FaceType::dirichlet ||
	    ownPrimary.faceTypes[face] == FaceType::analytical) {
      ltsSetup |= (1 << (face+4) );
    }
    // dynamic rupture faces are always global time stepping but operate on derivatives
    else if( ownPrimary.faceTypes[face] == FaceType::dynamicRupture ) {
      // face-neighbor provides GTS derivatives
      // face-neighbor provides derivatives
      ltsSetup |= ( 1 <<  face      );
      ltsSetup |= ( 1 << (face + 4) );

      // cell is required to provide derivatives for dynamic rupture
      ltsSetup |= ( 1 << 9 );

      if( copyLayer ) { // set the buffer invalid in copy layers
                     // TODO(unknown): Minor improvements possible: Non-DR MPI-neighbor for example
        ltsSetup |= ( 1 << 10 );
      }
    }
    // derive the LTS setup based on the cluster ids
    else {
      // neighboring cluster has a larger time step than this cluster
      if( ownSecondary.clusterId < secondary[face].get().clusterId ) {
        // neighbor delivers time derivatives
        ltsSetup |= ( 1 << face );

        // the cell-local buffer is used in LTS-fashion
        ltsSetup |= ( 1 << 10     );
      }
      // GTS relation
      else if( ownSecondary.clusterId == secondary[face].get().clusterId ) {
        ltsSetup |= ( 1 << (face + 4) );
      }

      // cell is required to provide derivatives
      if( ownSecondary.clusterId > secondary[face].get().clusterId ) {
        ltsSetup |= ( 1 << 9 );
      }
      // cell is required to provide a buffer
      else {
        ltsSetup |= ( 1 << 8 );
      }
    }

    // true lts buffer with gts required derivatives
    if( (ltsSetup >> 10)%2 == 1 && (ltsSetup >> 4)%16 != 0 ) {
      ltsSetup |= ( 1 << 9 );
    }
  }

  /*
   * Normalize for special case "free surface/dirichlet on derivatives":
   *   If a cell provides either buffers in a LTS fashion or derivatives only,
   *   the neighboring contribution of the boundary intergral is required to work on the cells derivatives.
   *   Free surface/dirichlet boundary conditions work on the cells DOFs in the neighboring contribution:
   *   Enable cell local derivatives in this case and mark that the "fake neighbor" provides derivatives.
   */
  for( unsigned int face = 0; face < 4; face++ ) {
    // check for special case free-surface/dirichlet requirements
    const bool isSpecialCase = ownPrimary.faceTypes[face] == FaceType::freeSurface ||
      ownPrimary.faceTypes[face] == FaceType::freeSurfaceGravity ||
      ownPrimary.faceTypes[face] == FaceType::dirichlet ||
      ownPrimary.faceTypes[face] == FaceType::analytical;
    if (isSpecialCase &&       // special case face
       ( (ltsSetup >> 10) % 2 == 1 ||             // lts fashion buffer
         (ltsSetup >> 8 ) % 2 == 0 )         ) {  // no buffer at all
      ltsSetup |= ( 1 << 9 );       // enable derivatives computation
      ltsSetup |= ( 1 << face);  // enable derivatives for the fake face neighbor
    }
  }

  return ltsSetup;
    }


/**
 * Normalizes the LTS setup for the special case "GTS on derivatives":
 *   If a face neighbor provides true buffers to cells with larger time steps,
 *   the local cell is required to operate on derivatives of this face neighbor.
 *
 *   Example:
 *        | own |  fn 1 |  fn 2 | fn 3 | fn 4 |
 *   local|  dt |    dt | 0.5dt |   dt |   dt |
 *   fn 4 |  dt | 0.5dt |   2dt |   dt |   dt |
 *         -----------------------------------
 *   In the example the local cell is connected via face 4 with to a GTS neighbor.
 *   Face neighbor 4 is required to deliver true buffers to its second face neighbor.
 *   It follows that the local cell has to operate on the derivatives of face neighbor 4.
 **/
    unsigned short correctLtsSetup(const CellLocalInformation& ownPrimary, const std::array<std::reference_wrapper<CellLocalInformation>, 4>& primary) {
        auto localLtsSetup = ownPrimary.ltsSetup;
// iterate over the face neighbors
  for( unsigned int face = 0; face < 4; face++ ) {
    // enforce derivatives if this is a "GTS on derivatives" relation
    if( (localLtsSetup >> (face + 4))%2 && (primary[face].get().ltsSetup >> 10)%2 == 1 ) {
      localLtsSetup |= (1 << face);
    }
  }
  return localLtsSetup;
    }

    unsigned short correctGhostLayer(int face, const CellLocalInformation& ownPrimary, const std::array<std::reference_wrapper<CellLocalInformation>, 4>& primary) {
        auto otherLtsSetup = primary[face].get().ltsSetup;
        
    }

    std::pair<std::array<std::reference_wrapper<const CellLocalInformation>, 4>, std::array<std::reference_wrapper<const SecondaryCellLocalInformation>, 4>> getNeighboringCellInformation(initializer::MemoryContainer& container, const CellLocalInformation& ownCellInformation, const SecondaryCellLocalInformation& ownSecondaryCellInformation) {
      std::vector<std::reference_wrapper<const CellLocalInformation>> cellInfo;
      std::vector<std::reference_wrapper<const SecondaryCellLocalInformation>> secondaryCellInfo;

      for (int face = 0; face < 4; ++face) {
          auto neighborConfig = ownCellInformation.neighborConfigIds[face];
          auto i2 = ownSecondaryCellInformation.faceNeighborIds[face];
          container.cluster.visitIdx(neighborConfig, [&](auto&& ltsview2) {
              auto* secondaryCellInformation2 = ltsview2.tree.var(ltsview2.lts.secondaryCellInformation);
              auto* cellInformation2 = ltsview2.tree.var(ltsview2.lts.cellInformation);
              cellInfo.emplace_back(std::cref(cellInformation2[i2]));
              secondaryCellInfo.emplace_back(std::cref(secondaryCellInformation2[i2]));
          });
      }
    }
} // namespace

namespace seissol::initializer::internal {
    void handleLtsSetup(MemoryContainer& container) {
      std::vector<unsigned short> copyLts(TODO);
container.cluster.visitLayers([&](auto&& layerview) {
        if (layerview.icg != Ghost) {
            auto* secondaryCellInformation = layerview.var(layerview.lts.secondaryCellInformation);
            auto* cellInformation = layerview.var(layerview.lts.cellInformation);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (int i = 0; i < layerview.layer.getNumberOfCells(); ++i) {
                std::array<std::reference_wrapper<const CellLocalInformation>, 4> cellInfo;
                std::array<std::reference_wrapper<const SecondaryCellLocalInformation>, 4> secondaryCellInfo;

                for (int face = 0; face < 4; ++face) {
                    auto neighborConfig = cellInformation[i].neighborConfigIds[face];
                    auto i2 = secondaryCellInformation[i].faceNeighborIds[face];
                    container.cluster.visitIdx(neighborConfig, [&](auto&& ltsview2) {
                        auto* secondaryCellInformation2 = ltsview2.var(ltsview2.lts.secondaryCellInformation);
                        auto* cellInformation2 = ltsview2.var(ltsview2.lts.cellInformation);
                        cellInfo[face] = std::cref(cellInformation2[i2]);
                        secondaryCellInfo[face] = std::cref(secondaryCellInformation2[i2]);
                    });
                }
                cellInformation[i].ltsSetup = determineLtsSetup(cellInformation[i], secondaryCellInformation[i], cellInfo, secondaryCellInfo, layerview.icg==Copy);
            }
        }
    });

    // do sync

    // pass 4: correct LTS setup
    container.cluster.visitLayers([&](auto&& layerview) {
        if (layerview.icg != Ghost) {
            auto* secondaryCellInformation = layerview.var(layerview.lts.secondaryCellInformation);
            auto* cellInformation = layerview.var(layerview.lts.cellInformation);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (int i = 0; i < layerview.layer.getNumberOfCells(); ++i) {
                std::array<std::reference_wrapper<const CellLocalInformation>, 4> cellInfo;
                std::array<std::reference_wrapper<const SecondaryCellLocalInformation>, 4> secondaryCellInfo;

                for (int face = 0; face < 4; ++face) {
                    auto neighborConfig = cellInformation[i].neighborConfigIds[face];
                    auto i2 = secondaryCellInformation[i].faceNeighborIds[face];
                    container.cluster.visitIdx(neighborConfig, [&](auto&& ltsview2) {
                        auto* secondaryCellInformation2 = ltsview2.var(ltsview2.lts.secondaryCellInformation);
                        auto* cellInformation2 = ltsview2.var(ltsview2.lts.cellInformation);
                        cellInfo[face] = std::cref(cellInformation2[i2]);
                        secondaryCellInfo[face] = std::cref(secondaryCellInformation2[i2]);
                    });
                }
                cellInformation[i].ltsSetup = correctLtsSetup(cellInformation[i], cellInfo);
            }
        }
    });

    // do sync
    container.cluster.visitLayers([&](auto&& layerview) {
        if (layerview.icg == Ghost) {
            auto* secondaryCellInformation = layerview.var(layerview.lts.secondaryCellInformation);
            auto* cellInformation = layerview.var(layerview.lts.cellInformation);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (int i = 0; i < layerview.layer.getNumberOfCells(); ++i) {
                std::array<std::reference_wrapper<const CellLocalInformation>, 4> cellInfo;
                std::array<std::reference_wrapper<const SecondaryCellLocalInformation>, 4> secondaryCellInfo;

                // TODO non-ghost face
                int face = 0; // TODO:
                auto neighborConfig = cellInformation[i].neighborConfigIds[face];
                auto i2 = secondaryCellInformation[i].faceNeighborIds[face];
                container.cluster.visitIdx(neighborConfig, [&](auto&& ltsview2) {
                    auto* secondaryCellInformation2 = ltsview2.var(ltsview2.lts.secondaryCellInformation);
                    auto* cellInformation2 = ltsview2.var(ltsview2.lts.cellInformation);
                    cellInfo[face] = std::cref(cellInformation2[i2]);
                    secondaryCellInfo[face] = std::cref(secondaryCellInformation2[i2]);
                });
                cellInformation[i].ltsSetup = correctGhostLayer(cellInformation[i], cellInfo);
            }
        }
    });
    }
} // namespace seissol::initializer::internal
