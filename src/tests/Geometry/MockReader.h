#pragma once

#include<array>
#include<Eigen/Dense>

#include "Geometry/MeshReader.h"

namespace seissol
{
  class MockReader : public MeshReader
  {
    public:
      MockReader(std::array<Eigen::Vector3d, 4> vertices) : MeshReader(0){
        m_vertices.resize(4);
        for (int i = 0; i < 4; i++) {
          std::copy(vertices[i].data(), vertices[i].data()+3, m_vertices.at(i).coords);
          m_vertices.at(i).elements = {1};
        }


        m_elements.resize(1);
        m_elements.at(0).vertices[0] = 0;
        m_elements.at(0).vertices[1] = 1;
        m_elements.at(0).vertices[2] = 2;
        m_elements.at(0).vertices[3] = 3;

      }
  };
}
