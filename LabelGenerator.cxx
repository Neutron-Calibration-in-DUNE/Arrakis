/**
 * @file LabelGenerator.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-07-21
 */
#include "LabelGenerator.h"

namespace arrakis
{
    LabelGenerator::LabelGenerator()
    {
        mLabelTree = mTFileService->make<TTree>("labels", "labels");
        mLabelTree->Branch("u1_gamma_id", &mLabels.u1_gamma_id);
        mLabelTree->Branch("u1_gamma_type", &mLabels.u1_gamma_type);
        mLabelTree->Branch("u1_ancestor_pdg", &mLabels.u1_ancestor_pdg);
        mLabelTree->Branch("u1_ancestor_id", &mLabels.u1_ancestor_id);
        mLabelTree->Branch("v1_gamma_id", &mLabels.v1_gamma_id);
        mLabelTree->Branch("v1_gamma_type", &mLabels.v1_gamma_type);
        mLabelTree->Branch("v1_ancestor_pdg", &mLabels.v1_ancestor_pdg);
        mLabelTree->Branch("v1_ancestor_id", &mLabels.v1_ancestor_id);
        mLabelTree->Branch("z1_gamma_id", &mLabels.z1_gamma_id);
        mLabelTree->Branch("z1_gamma_type", &mLabels.z1_gamma_type);
        mLabelTree->Branch("z1_ancestor_pdg", &mLabels.z1_ancestor_pdg);
        mLabelTree->Branch("z1_ancestor_id", &mLabels.z1_ancestor_id);

        mLabelTree->Branch("u2_gamma_id", &mLabels.u2_gamma_id);
        mLabelTree->Branch("u2_gamma_type", &mLabels.u2_gamma_type);
        mLabelTree->Branch("v2_gamma_id", &mLabels.v2_gamma_id);
        mLabelTree->Branch("v2_gamma_type", &mLabels.v2_gamma_type);
        mLabelTree->Branch("z2_gamma_id", &mLabels.z2_gamma_id);
        mLabelTree->Branch("z2_gamma_type", &mLabels.z2_gamma_type);
    }

    LabelGenerator::~LabelGenerator()
    {}

    void LabelGenerator::ResetLabels()
    {
        mLabels.u1_gamma_id.clear();
        mLabels.u1_gamma_type.clear();
        mLabels.u1_ancestor_pdg.clear();
        mLabels.u1_ancestor_id.clear();

        mLabels.v1_gamma_id.clear();
        mLabels.v1_gamma_type.clear();
        mLabels.v1_ancestor_pdg.clear();
        mLabels.v1_ancestor_id.clear();

        mLabels.z1_gamma_id.clear();
        mLabels.z1_gamma_type.clear();
        mLabels.z1_ancestor_pdg.clear();
        mLabels.z1_ancestor_id.clear();

        mLabels.u2_gamma_id.clear();
        mLabels.u2_gamma_type.clear();
        mLabels.v2_gamma_id.clear();
        mLabels.v2_gamma_type.clear();
        mLabels.z2_gamma_id.clear();
        mLabels.z2_gamma_type.clear();
    }

    int LabelGenerator::getLargestEnergyTrackID(std::vector<int> trackID_vec, std::vector<double> energy_vec)
    {
        int track_id = 0;
        double track_energy;
        int track_index = 0;
        track_energy = energy_vec[0];
        for (size_t i = 0; i < energy_vec.size(); i++)
        {
            if (energy_vec[i] > track_energy)
            {
                track_energy = energy_vec[i];
                track_index = i;
            }
        }
        track_id = trackID_vec[track_index];
        return track_id;
    }

    void LabelGenerator::processEvent(
        ParticleTree particleTree,
        GammaTable gammaTable,
        ArrayGenerator arrayGenerator
    )
    {
        ResetLabels();

        EventArray eventArray = arrayGenerator.getArray();
        std::vector<Gamma> gammas = gammaTable.getGammaTable();
        std::map<Int_t, Int_t> gammaTableIndex = gammaTable.getGammaTableIndex();
        std::vector<Int_t> gammaTableTrackIDs = gammaTable.getGammaTableTrackIDs();

        for(int j = 0; j< (int) eventArray.u1_tdc.size(); j++)
        {
            if(eventArray.u1_track_ids[j].empty() == 0 && eventArray.u1_energy[j].empty() == 0)
            {
                int temp_trackID = getLargestEnergyTrackID(eventArray.u1_track_ids[j], eventArray.u1_energy[j]);

                // check if track id is in the gamma table
                bool track_id_exists = false;
                Int_t gammaIndex = 0;
                for (size_t i = 0; i < gammaTableTrackIDs.size(); i++)
                {
                    if (gammaTableTrackIDs[i] == temp_trackID) 
                    {
                        track_id_exists = true;
                        gammaIndex = gammaTableIndex[temp_trackID];
                    }
                }
                if (track_id_exists)
                {
                    mLabels.u1_gamma_id.emplace_back(gammas[gammaIndex].track_id);
                    mLabels.u1_gamma_type.emplace_back(gammas[gammaIndex].energy);
                }
                else
                {
                    mLabels.u1_gamma_id.emplace_back(-1);
                    mLabels.u1_gamma_type.emplace_back(0);
                }
                mLabels.u1_ancestor_pdg.emplace_back(particleTree.GetAncestorPDG(temp_trackID));
                mLabels.u1_ancestor_id.emplace_back(particleTree.GetAncestorTrackID(temp_trackID));
            }
            else
            {
                mLabels.u1_gamma_id.emplace_back(-2);
                mLabels.u1_gamma_type.emplace_back(-999);
                mLabels.u1_ancestor_pdg.emplace_back(0);
                mLabels.u1_ancestor_id.emplace_back(0);
            }
		}
        mLabelTree->Fill();
    }
}