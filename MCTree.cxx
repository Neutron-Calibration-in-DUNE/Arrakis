/**
 * @file MCTree.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-02-20
 */
#include "MCTree.h"

namespace arrakis
{
    MCTree::MCTree()
    {
    }

    MCTree::~MCTree()
    {
    }

    void MCTree::ProcessEventMC(
        ParticleMaps* particle_maps,
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
        const art::ValidHandle<std::vector<sim::SimEnergyDeposit>>& mcEnergyDeposits
    )
    {
        
    }
}