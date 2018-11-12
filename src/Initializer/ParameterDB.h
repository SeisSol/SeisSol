/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2017, SeisSol Group
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
 * Setup of SeisSol's cell local matrices.
 **/

#ifndef INITIALIZER_PARAMETERDB_H_
#define INITIALIZER_PARAMETERDB_H_

#include <string>
#include <unordered_map>
#include <Geometry/MeshReader.h>
#include <easi/Query.h>

#ifndef PUML_PUML_H
namespace PUML
{
	class TETPUML;
}
#endif // PUML_PUML_H

namespace seissol {
  namespace initializers {
    class QueryGenerator;
    class ElementBarycentreGenerator;
    class ElementBarycentreGeneratorPUML;
    class FaultBarycentreGenerator;
    class FaultGPGenerator;
    class ParameterDB;
  }
}

namespace easi {
  class Component;
};

class seissol::initializers::QueryGenerator {
public:
  virtual easi::Query generate() const = 0;
};

class seissol::initializers::ElementBarycentreGenerator : public seissol::initializers::QueryGenerator {
public:
  ElementBarycentreGenerator(MeshReader const& meshReader) : m_meshReader(meshReader) {}
  virtual easi::Query generate() const;
private:
  MeshReader const& m_meshReader;
};

#ifdef USE_HDF
class seissol::initializers::ElementBarycentreGeneratorPUML : public seissol::initializers::QueryGenerator {
public:
  ElementBarycentreGeneratorPUML(PUML::TETPUML const& mesh) : m_mesh(mesh) {}
  virtual easi::Query generate() const;
private:
  PUML::TETPUML const& m_mesh;
};
#endif

class seissol::initializers::FaultBarycentreGenerator : public seissol::initializers::QueryGenerator {
public:
  FaultBarycentreGenerator(MeshReader const& meshReader, unsigned numberOfPoints) : m_meshReader(meshReader), m_numberOfPoints(numberOfPoints) {}
  virtual easi::Query generate() const;

private:
  MeshReader const& m_meshReader;
  unsigned m_numberOfPoints;
};

class seissol::initializers::FaultGPGenerator : public seissol::initializers::QueryGenerator {
public:
  FaultGPGenerator(MeshReader const& meshReader, double (*points)[2], unsigned numberOfPoints) : m_meshReader(meshReader), m_points(points), m_numberOfPoints(numberOfPoints) {}
  virtual easi::Query generate() const;
private:
  MeshReader const& m_meshReader;
  double (*m_points)[2];
  unsigned m_numberOfPoints;
};

class seissol::initializers::ParameterDB {
public:  
  void addParameter(std::string const& parameter, double* memory, unsigned stride = 1) { m_parameters[parameter] = std::make_pair(memory, stride); }
  void evaluateModel(std::string const& fileName, QueryGenerator const& queryGen);
  static bool faultParameterizedByTraction(std::string const& fileName);
  
private:
  static easi::Component* loadModel(std::string const& fileName);
  std::unordered_map<std::string, std::pair<double*, unsigned>> m_parameters;
};

#endif
