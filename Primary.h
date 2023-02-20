/**
 * @file Primary.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-02-20
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
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft includes
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCGeneratorInfo.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"

// necessary ROOT libraries
#include <TTree.h>

#include "Generators.h"
#include "Logger.h"
#include "Node.h"

namespace arrakis
{
    class Primary : public Node
    {
    public:
        Primary();
        ~Primary();
        
        Primary(GeneratorLabel label, simb::MCParticle& particle);

    private:
        
    };
}