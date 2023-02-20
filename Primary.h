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
#include "TrajectoryStep.h"

namespace arrakis
{
    class Primary : public Node
    {
    public:
        Primary();
        ~Primary();
        
        Primary(GeneratorLabel label, simb::MCParticle& particle);

        GeneratorLabel GetGeneratorLabel()  { return mGeneratorLabel; }
        Int_t TrackID() { return mTrackID; }
        Int_t PdgCode() { return mPdgCode; }
        std::string InitProcess()   { return mInitProcess; }
        Double_t InitEnergy()       { return mInitEnergy; }
        Double_t InitX()            { return mInitX; }
        Double_t InitY()            { return mInitY; }
        Double_t InitZ()            { return mInitZ; }
        std::string EndProcess()    { return mEndProcess; }
        Double_t EndEnergy()        { return mEndEnergy; }
        Double_t EndX()             { return mEndX; }
        Double_t EndY()             { return mEndY; }
        Double_t EndZ()             { return mEndZ; }

    private:
        GeneratorLabel mGeneratorLabel = {kNone};
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

        // First trajectory point
        TrajectoryStep* mTrajectoryStep = {0};
    };
}