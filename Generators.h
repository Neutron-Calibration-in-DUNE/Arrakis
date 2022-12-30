/**
 * @file Generators.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-12/21
 */
#pragma once
#include <vector>
#include <string>
#include <algorithm>
#include <map>

// special utility includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCGeneratorInfo.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCTrajectory.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/Simulation/sim.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/RawDigit.h"

// necessary ROOT libraries
#include <TTree.h>

namespace arrakis
{
    struct Generator
    {
        art::InputTag input_tag;
        std::string label;
        art::ValidHandle<std::vector<simb::MCTruth>> truth;

        Generator(
            art::InputTag _input_tag,
            std::string _label, 
            art::ValidHandle<std::vector<simb::MCTruth>> _truth
        )
        : input_tag(_input_tag)
        , label(_label)
        , truth(_truth)
        {
        }
    };

    class Generators
    {
    public:
        Generators();
        ~Generators();

        void ProcessEvent(
            std::vector<art::InputTag> inputTags,
            std::vector<std::string> labels,
            std::vector<art::ValidHandle<std::vector<simb::MCTruth>>> mcTruth
        );

        const std::vector<Generator>& GetGenerators() { return mGenerators; }

    private:
        std::vector<Generator> mGenerators;
    };
}