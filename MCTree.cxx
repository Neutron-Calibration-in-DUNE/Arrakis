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

    void MCTree::ResetEvent()
    {
        mPrimaries.clear();
    }

    void MCTree::ProcessEventMC(
        ParticleMaps* particle_maps,
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
        const art::ValidHandle<std::vector<sim::SimEnergyDeposit>>& mcEnergyDeposits
    )
    {
        if (!mcParticles.isValid()) {
            Logger::GetInstance("mc_tree")->error("MCParticles handle is not valid!");
            return;
        }
        if (!mcEnergyDeposits.isValid()) {
            Logger::GetInstance("mc_tree")->error("SimEnergyDeposits handle is not valid!");
            return;
        }
        ResetEvent();
        Logger::GetInstance("mc_tree")->trace("processing " + std::to_string((*mcParticles).size()) + " MCParticles");
        // Loop through every particle and save
        // information to the primary that was
        // generated.
        for (auto particle : *mcParticles)
        {
            // If the particle is a primary, make
            // a new entry in mPrimaries.
            if(particle.Mother() == 0) 
            {
                mPrimaries.emplace_back(new Primary(
                    particle_maps->GetGeneratorLabel(particle.TrackId()),
                    particle
                ));
            }
            // Otherwise, find the associated primary
            // using the ancestor map.
            else
            {
                // Int_t primary_index = FindPrimary(
                //     particle_maps->GetAncestorTrackID(particle.TrackId())
                // );
                // if(primary_index == -1) 
                // {
                //     Logger::GetInstance("mc_tree")->warning(
                //         "could not find primary with track id " + 
                //         std::to_string(particle_maps->GetAncestorTrackID(particle.TrackId())) + 
                //         " from particle ancestor with track id " + 
                //         std::to_string(particle.TrackId())
                //     );
                //     continue;
                // }
                // mPrimaries[primary_index].AddDaughter(
                //     particle, particle_maps->GetAncestorLevel(particle.TrackId())
                // );
            }
        }
    }
}