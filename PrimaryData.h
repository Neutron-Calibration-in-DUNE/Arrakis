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
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// necessary ROOT libraries
#include <TTree.h>

#include "ParticleMaps.h"

namespace arrakis
{
    struct Primary
    {
        Int_t track_id = {0};
        Int_t pdg = {0};

        std::string init_process = {""};
        Double_t init_energy = {0};
        Double_t init_x = {0};
        Double_t init_y = {0};
        Double_t init_z = {0};
        
        std::string end_process = {""};
        Double_t end_energy = {0};
        Double_t end_x = {0};
        Double_t end_y = {0};
        Double_t end_z = {0};

        std::vector<Int_t> daughter_ids = {};
        std::vector<Int_t> daughter_level = {};

        std::vector<std::string> daughter_init_process = {};
        std::vector<Double_t> daughter_init_energy = {};
        std::vector<Double_t> daughter_init_x = {};
        std::vector<Double_t> daughter_init_y = {};
        std::vector<Double_t> daughter_init_z = {};

        std::vector<std::string> daughter_end_process = {};
        std::vector<Double_t> daughter_end_energy = {};
        std::vector<Double_t> daughter_end_x = {};
        std::vector<Double_t> daughter_end_y = {};
        std::vector<Double_t> daughter_end_z = {};

        // std::vector<std::vector<Double_t>> daughter_edep_energy = {};
        // std::vector<std::vector<Double_t>> daughter_edep_x = {};
        // std::vector<std::vector<Double_t>> daughter_edep_y = {};
        // std::vector<std::vector<Double_t>> daughter_edep_z = {};
        // std::vector<std::vector<Int_t>> daughter_edep_num_electrons = {};
        // std::vector<std::vector<Int_t>> daughter_edep_num_photons = {};
        // Int_t num_edep_points = 0;
    };

    class PrimaryData
    {
    public:
        PrimaryData();
        ~PrimaryData();

        void ResetEvent();
        void ProcessEvent(
            ParticleMaps particle_maps,
            const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
            const art::ValidHandle<std::vector<simb::EnergyDeposit>>& mcEnergyDeposits
        );

        Int_t FindPrimary(Int_t track_id);

    private:
        art::ServiceHandle<art::TFileService> mTFileService;
        TTree *mTTree;

        std::vector<Primary> mPrimaries;
        Primary mPrimary;

    };
}