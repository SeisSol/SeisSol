/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Vishal Sontakke (vishal.sontakke AT tum.de)
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
 */

#ifndef POST_PROCESSOR_H
#define POST_PROCESSOR_H

#include <vector>
#include <Initializer/typedefs.hpp>
#include <Initializer/tree/Layer.hpp>
#include <Initializer/tree/LTSTree.hpp>
#include <Initializer/preProcessorMacros.fpp>

namespace seissol
{

namespace writer
{

class PostProcessor {
private:
    bool m_integrationMask[9];
    int m_numberOfVariables;
    std::vector<int> m_integerMap;
    seissol::initializers::Variable<real> m_integrals;
public:
    PostProcessor (): m_numberOfVariables(0), m_integerMap(0L) {
        for (size_t i = 0; i < 9; i++) {
            m_integrationMask[i] = false;
        }
    }
    virtual ~PostProcessor () {}
    void integrateQuantities(const double i_timestep,
    	seissol::initializers::Layer& i_layerData, const unsigned int l_cell,
    	const double * const i_dofs);
    void setIntegrationMask(const int * const i_integrationMask);
    int getNumberOfVariables();
    void getIntegrationMask(bool* transferTo);
    void allocateMemory(seissol::initializers::LTSTree* ltsTree);
    const real* getIntegrals(seissol::initializers::LTSTree* ltsTree);
};

}

}

#endif // POST_PROCESSOR_H
