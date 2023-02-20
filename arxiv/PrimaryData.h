/**
 * @file PrimaryData.h
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
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/Simulation/sim.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/RawDigit.h"

// necessary ROOT libraries
#include <TTree.h>

#include "DetectorGeometry.h"
#include "Logger.h"
#include "ParticleMaps.h"
#include "Primary.h"
#include "Junk.h"

namespace arrakis
{
    /**
     * @brief This class holds a large amount of MC truth data
     * for each primary particle that is generated in Geant4.
     * The truth data includes MCParticle, SimEnergyDeposit,
     * SimChannel and RawDigit information for the primary
     * and all of its daughters.
     * 
     * The data can be saved to a TTree, although this can 
     * be quite a large amount of information.  
     * 
     */
    class PrimaryData
    {
    public:
        PrimaryData(
            bool SavePrimaryData, bool SavePrimaryDataEdeps,
            bool SavePrimaryDataRawTPC,
            Double_t ADCThreshold
        );
        ~PrimaryData();

        void ResetEvent();
        void ProcessEventMC(
            ParticleMaps* particle_maps,
            const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
            const art::ValidHandle<std::vector<sim::SimEnergyDeposit>>& mcEnergyDeposits
        );
        void ProcessEventDetectorSimulation(
            ParticleMaps* particle_maps,
            detinfo::DetectorClocksData const& clockData,
            const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
            const art::ValidHandle<std::vector<raw::RawDigit>>& rawTPC
        );
        void FillTTree();
        Int_t FindPrimary(Int_t track_id);
        void FindDetectorProcess(
            detinfo::DetectorClocksData const& clockData,
            Int_t primary_index, Int_t track_id, 
            Double_t energy, unsigned int tdc
        );

        void PrintPrimaryEnergyDepositions();
        void PrintDaughterEnergyDepositions();

        std::vector<Primary> GetPrimaries() { return mPrimaries; }
        Junk GetJunk() { return mJunk; }

    private:
        bool mSavePrimaryData = {false};
        bool mSavePrimaryDataEdeps = {false};
        bool mSavePrimaryDataRawTPC = {false};

        Double_t mADCThreshold = {4.0};

        art::ServiceHandle<art::TFileService> mTFileService;
        TTree *mTTree;

        std::vector<Primary> mPrimaries;
        Junk mJunk;
        Primary mPrimary;

    };
}