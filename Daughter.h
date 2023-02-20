/**
 * @file Daughter.h
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
#include "TrajectoryStep.h"

namespace arrakis
{
    class Daughter : public Node
    {
    public:.
        Daughter();
        ~Daughter();
        
        Daughter(simb::MCParticle& particle);

        const Int_t& TrackID() { return mTrackID; }
        const Int_t& PdgCode() { return mPdgCode; }
        const std::string& InitProcess() const { return mInitProcess; }
        const Double_t& InitEnergy() const     { return mInitEnergy; }
        const Double_t& InitX() const          { return mInitX; }
        const Double_t& InitY() const          { return mInitY; }
        const Double_t& InitZ() const          { return mInitZ; }
        const std::string& EndProcess() const  { return mEndProcess; }
        const Double_t& EndEnergy() const      { return mEndEnergy; }
        const Double_t& EndX() const           { return mEndX; }
        const Double_t& EndY() const           { return mEndY; }
        const Double_t& EndZ() const           { return mEndZ; }
    private:
        Int_t mTrackID = {0};
        Int_t mPdgCode = {0};

        // MC Particle info.
        std::string mInitProcess = {""};
        Double_t mInitEnergy = {0};
        Double_t mInitT = {0};
        Double_t mInitX = {0};
        Double_t mInitY = {0};
        Double_t mInitZ = {0};
        
        std::string mEndProcess = {""};
        Double_t mEndEnergy = {0};
        Double_t mEndT = {0};
        Double_t mEndX = {0};
        Double_t mEndY = {0};
        Double_t mEndZ = {0};

    };
}