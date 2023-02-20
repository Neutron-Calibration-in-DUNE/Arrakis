/**
 * @file TrajectoryStep.h
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

#include "Daughter.h"
#include "DetectorGeometry.h"
#include "EnergyDeposition.h"
#include "Generators.h"
#include "Logger.h"
#include "Node.h"

namespace arrakis
{
    class TrajectoryStep : public Node
    {
    public:
        TrajectoryStep();
        ~TrajectoryStep();
        
        TrajectoryStep(
            size_t step, simb::MCParticle& particle, 
            simb::MCTrajectory& trajectory, 
            std::vector<std::string>& trajectory_processes
        );

        const std::string& Process() const   { return mProcess; }
        const DetectorVolume& GetDetectorVolume() const  { return mDetectorVolume; }
        const Double_t& T() const { return mT; }
        const Double_t& X() const { return mX; }
        const Double_t& Y() const { return mY; }
        const Double_t& Z() const { return mZ; }

        Int_t NumberOfDaughters()  { return mDaughters.size(); }

    private:
        std::string mProcess = {"not_defined"};
        DetectorVolume mDetectorVolume;
        Double_t mT = {0};
        Double_t mX = {0};
        Double_t mY = {0};
        Double_t mZ = {0};

        TrajectoryStep* mNextTrajectoryStep = {0};
        EnergyDeposition* mEnergyDeposition = {0};
        std::vector<Daughter*> mDaughters = {};
    };
}