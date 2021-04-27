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

#include "PostProcessor.h"
#include "SeisSol.h"

void seissol::writer::PostProcessor::integrateQuantities(const double i_timestep,
	seissol::initializers::Layer& i_layerData, const unsigned int l_cell,
	const double * const i_dofs) {

	real *integrals = i_layerData.var(m_integrals);
	for (int i = 0; i < m_numberOfVariables; i++) {
		integrals[l_cell*m_numberOfVariables+i] += i_dofs[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*m_integerMap[i]]*i_timestep;
	}
}

void seissol::writer::PostProcessor::setIntegrationMask(const int * const i_integrationMask) {
	unsigned int nextId = 0;
	for (int i = 0; i < 9; i++) {
		m_integrationMask[i] = (i_integrationMask[i] > 0);
		if (m_integrationMask[i]) {
			m_integerMap.push_back(i);
			m_numberOfVariables++;
			nextId++;
		}
	}
	m_integrals.count = m_numberOfVariables;
}

int seissol::writer::PostProcessor::getNumberOfVariables() {
	return m_numberOfVariables;
}

void seissol::writer::PostProcessor::getIntegrationMask(bool* transferTo) {
	for(int i = 0; i < 9; i++) {
		transferTo[i] = m_integrationMask[i];
	}
}

void seissol::writer::PostProcessor::allocateMemory(seissol::initializers::LTSTree* ltsTree) {
	ltsTree->addVar( m_integrals, seissol::initializers::LayerMask(Ghost), PAGESIZE_HEAP,
      seissol::memory::Standard );
}

const real* seissol::writer::PostProcessor::getIntegrals(seissol::initializers::LTSTree* ltsTree) {
	if (m_numberOfVariables == 0) {
		return 0L;
	} else {
		return ltsTree->var(m_integrals);
	}
}
