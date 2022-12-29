/**
 * @file SoloPointCloudGenerator.h
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
#include "lardataobj/RawData/RawDigit.h"

// necessary ROOT libraries
#include <TTree.h>

#include "ParticleMaps.h"
#include "PrimaryData.h"

namespace arrakis
{
    enum SoloPointCloudLabel
    {
        "None" = 0,
        "NeutronElastic" = 1,
        "NeutronCapture" = 2,
        "NeutronGamma_4.75" = 3,
        "NeutronGamma_1.18" = 4
    };

    struct SoloPointCloud
    {
        Int_t event_id = {-1};
        Int_t point_cloud_id = {-1};
        std::vector<Double_t> channel = {};
        std::vector<Double_t> tdc = {};
        std::vector<Double_t> adc = {};
        std::vector<Double_t> energy = {};
        enum SoloPointCloudLabel label = {};

        bool all_deposited = false;
        bool all_lar = false;
        bool same_apa = false;
        Double_t lar_edep_fraction = {0};
    };


    class SoloPointCloudGenerator
    {
    public:
        SoloPointCloudGenerator();
        ~SoloPointCloudGenerator();

        void ProcessEvent(ParticleMaps particle_maps, PrimaryData primary_data);
    private:
        art::ServiceHandle<art::TFileService> mTFileService;
        TTree *mTTree;

        Double_t mEdepEnergyThreshold = {0.01};

        std::vector<SoloPointCloud> mSoloPointClouds;
        SoloPointCloud mSoloPointCloud;
    };
}