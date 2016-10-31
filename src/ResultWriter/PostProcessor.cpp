#include "PostProcessor.h"
#include "SeisSol.h"

void seissol::writer::PostProcessor::integrateQuantities(const double i_timestep,
	seissol::initializers::Layer& i_layerData, const unsigned int l_cell,
	const double * const i_dofs) {

	real (*integrals) = i_layerData.var(m_integrals);
	for (unsigned int i = 0; i < m_numberOfVariables; i++) {
		integrals[l_cell*m_numberOfVariables+i] += i_dofs[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*m_integerMap[i]]*i_timestep;
	}
}

void seissol::writer::PostProcessor::setIntegrationMask(const int * const i_integrationMask) {
	unsigned int nextId = 0;
	for (unsigned int i = 0; i < 9; i++) {
		m_integrationMask[i] = (bool)i_integrationMask[i];
		if (m_integrationMask[i]) {
			m_integerMap.push_back(nextId);
			m_numberOfVariables++;
			nextId++;
		}
	}
	m_integrals.count = m_numberOfVariables;
}

int seissol::writer::PostProcessor::getNumberOfVariables() {
	return m_numberOfVariables;
}

bool* seissol::writer::PostProcessor::getIntegrationMask() {
	return m_integrationMask;
}

void seissol::writer::PostProcessor::allocateMemory(seissol::initializers::LTSTree* ltsTree) {
	ltsTree->addVar( m_integrals, seissol::initializers::LayerMask(Ghost), PAGESIZE_HEAP,
      seissol::memory::Standard );
}

const double* seissol::writer::PostProcessor::getIntegrals(seissol::initializers::LTSTree* ltsTree) {
	return reinterpret_cast<const double*>(ltsTree->var(m_integrals));
}
