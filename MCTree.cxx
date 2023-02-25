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
    namespace mctree
    {
        MCTree::MCTree()
        {
        }

        MCTree::~MCTree()
        {
        }

        void MCTree::ResetEvent()
        {
            // if(mPrimaries.size() == 0) { 
            //     return;
            // }
            // mPrimaries.clear();
        }

        void MCTree::ProcessEvent(const Parameters& config, art::Event const& event)
        {
            ProcessMCParticles(config, event);
            // if (!mcParticles.isValid()) {
            //     Logger::GetInstance("mc_tree")->error("MCParticles handle is not valid!");
            //     return;
            // }
            // if (!mcEnergyDeposits.isValid()) {
            //     Logger::GetInstance("mc_tree")->error("SimEnergyDeposits handle is not valid!");
            //     return;
            // }
            // //ResetEvent();
            // Logger::GetInstance("mc_tree")->trace("processing " + std::to_string((*mcParticles).size()) + " MCParticles");
            // // Loop through every particle and save
            // // information to the primary that was
            // // generated.
            // for (auto particle : *mcParticles)
            // {
            //     // If the particle is a primary, make
            //     // a new entry in mPrimaries.
            //     if(particle.Mother() == 0) 
            //     {
            //         mPrimaries.emplace_back(new Primary(
            //             particle_maps->GetGeneratorLabel(particle.TrackId()),
            //             particle
            //         ));
            //     }
            //     // Otherwise, find the associated primary
            //     // using the ancestor map.
            //     else
            //     {
            //         // Int_t primary_index = FindPrimary(
            //         //     particle_maps->GetAncestorTrackID(particle.TrackId())
            //         // );
            //         // if(primary_index == -1) 
            //         // {
            //         //     Logger::GetInstance("mc_tree")->warning(
            //         //         "could not find primary with track id " + 
            //         //         std::to_string(particle_maps->GetAncestorTrackID(particle.TrackId())) + 
            //         //         " from particle ancestor with track id " + 
            //         //         std::to_string(particle.TrackId())
            //         //     );
            //         //     continue;
            //         // }
            //         // mPrimaries[primary_index].AddDaughter(
            //         //     particle, particle_maps->GetAncestorLevel(particle.TrackId())
            //         // );
            //     }
            // }
        }
        void MCTree::ProcessMCParticles(const Parameters& config, art::Event const& event)
        {
            Logger::GetInstance("mctree")->trace(
                "Creating primary nodes..."
            );
            auto mc_data = mcdata::MCData::GetInstance();
            auto mc_particles = *mc_data->GetMCParticles();
            std::cout << mc_particles.size() << std::endl;
            for(size_t ii = 0; ii < mc_particles.size(); ii++)
            {
                std::cout << mc_particles[ii].TrackId() << std::endl;
                // If the particle is a primary, make
                // a new entry in mPrimaries.
                if(mc_particles[ii].Mother() == 0) 
                {
                    std::shared_ptr<Node> primary = CreatePrimary(mc_particles[ii], ii);
                    sPrimaries[mc_particles[ii].TrackId()] = primary;
                }
            }
        }
        std::shared_ptr<Node> MCTree::CreatePrimary(const simb::MCParticle& particle, Int_t index)
        {
            std::shared_ptr<Node> primary = std::make_shared<Node>(PrimaryData(index));

            return primary;
        }
    }
}