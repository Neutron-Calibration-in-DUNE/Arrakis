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
        mLabelTree->Branch("u2_ancestor_pdg", &mLabels.u2_ancestor_pdg);
        mLabelTree->Branch("u2_ancestor_id", &mLabels.u2_ancestor_id);

        mLabelTree->Branch("v2_gamma_id", &mLabels.v2_gamma_id);
        mLabelTree->Branch("v2_gamma_type", &mLabels.v2_gamma_type);
        mLabelTree->Branch("v2_ancestor_pdg", &mLabels.v2_ancestor_pdg);
        mLabelTree->Branch("v2_ancestor_id", &mLabels.v2_ancestor_id);

        mLabelTree->Branch("z2_gamma_id", &mLabels.z2_gamma_id);
        mLabelTree->Branch("z2_gamma_type", &mLabels.z2_gamma_type);
        mLabelTree->Branch("z2_ancestor_pdg", &mLabels.z2_ancestor_pdg);
        mLabelTree->Branch("z2_ancestor_id", &mLabels.z2_ancestor_id);
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
        mLabels.u2_ancestor_pdg.clear();
        mLabels.u2_ancestor_id.clear();

        mLabels.v2_gamma_id.clear();
        mLabels.v2_gamma_type.clear();
        mLabels.v2_ancestor_pdg.clear();
        mLabels.v2_ancestor_id.clear();

        mLabels.z2_gamma_id.clear();
        mLabels.z2_gamma_type.clear();
        mLabels.z2_ancestor_pdg.clear();
        mLabels.z2_ancestor_id.clear();
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
        for(int j = 0; j< (int) eventArray.v1_tdc.size(); j++)
        {
            if(eventArray.v1_track_ids[j].empty() == 0 && eventArray.v1_energy[j].empty() == 0)
            {
                int temp_trackID = getLargestEnergyTrackID(eventArray.v1_track_ids[j], eventArray.v1_energy[j]);

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
                    mLabels.v1_gamma_id.emplace_back(gammas[gammaIndex].track_id);
                    mLabels.v1_gamma_type.emplace_back(gammas[gammaIndex].energy);
                }
                else
                {
                    mLabels.v1_gamma_id.emplace_back(-1);
                    mLabels.v1_gamma_type.emplace_back(0);
                }
                mLabels.v1_ancestor_pdg.emplace_back(particleTree.GetAncestorPDG(temp_trackID));
                mLabels.v1_ancestor_id.emplace_back(particleTree.GetAncestorTrackID(temp_trackID));
            }
            else
            {
                mLabels.v1_gamma_id.emplace_back(-2);
                mLabels.v1_gamma_type.emplace_back(-999);
                mLabels.v1_ancestor_pdg.emplace_back(0);
                mLabels.v1_ancestor_id.emplace_back(0);
            }
		}
        for(int j = 0; j< (int) eventArray.z1_tdc.size(); j++)
        {
            if(eventArray.z1_track_ids[j].empty() == 0 && eventArray.z1_energy[j].empty() == 0)
            {
                int temp_trackID = getLargestEnergyTrackID(eventArray.z1_track_ids[j], eventArray.z1_energy[j]);

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
                    mLabels.z1_gamma_id.emplace_back(gammas[gammaIndex].track_id);
                    mLabels.z1_gamma_type.emplace_back(gammas[gammaIndex].energy);
                }
                else
                {
                    mLabels.z1_gamma_id.emplace_back(-1);
                    mLabels.z1_gamma_type.emplace_back(0);
                }
                mLabels.z1_ancestor_pdg.emplace_back(particleTree.GetAncestorPDG(temp_trackID));
                mLabels.z1_ancestor_id.emplace_back(particleTree.GetAncestorTrackID(temp_trackID));
            }
            else
            {
                mLabels.z1_gamma_id.emplace_back(-2);
                mLabels.z1_gamma_type.emplace_back(-999);
                mLabels.z1_ancestor_pdg.emplace_back(0);
                mLabels.z1_ancestor_id.emplace_back(0);
            }
		}
        for(int j = 0; j< (int) eventArray.u2_tdc.size(); j++)
        {
            if(eventArray.u2_track_ids[j].empty() == 0 && eventArray.u2_energy[j].empty() == 0)
            {
                int temp_trackID = getLargestEnergyTrackID(eventArray.u2_track_ids[j], eventArray.u2_energy[j]);

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
                    mLabels.u2_gamma_id.emplace_back(gammas[gammaIndex].track_id);
                    mLabels.u2_gamma_type.emplace_back(gammas[gammaIndex].energy);
                }
                else
                {
                    mLabels.u2_gamma_id.emplace_back(-1);
                    mLabels.u2_gamma_type.emplace_back(0);
                }
                mLabels.u2_ancestor_pdg.emplace_back(particleTree.GetAncestorPDG(temp_trackID));
                mLabels.u2_ancestor_id.emplace_back(particleTree.GetAncestorTrackID(temp_trackID));
            }
            else
            {
                mLabels.u2_gamma_id.emplace_back(-2);
                mLabels.u2_gamma_type.emplace_back(-999);
                mLabels.u2_ancestor_pdg.emplace_back(0);
                mLabels.u2_ancestor_id.emplace_back(0);
            }
		}
        for(int j = 0; j< (int) eventArray.v2_tdc.size(); j++)
        {
            if(eventArray.v2_track_ids[j].empty() == 0 && eventArray.v2_energy[j].empty() == 0)
            {
                int temp_trackID = getLargestEnergyTrackID(eventArray.v2_track_ids[j], eventArray.v2_energy[j]);

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
                    mLabels.v2_gamma_id.emplace_back(gammas[gammaIndex].track_id);
                    mLabels.v2_gamma_type.emplace_back(gammas[gammaIndex].energy);
                }
                else
                {
                    mLabels.v2_gamma_id.emplace_back(-1);
                    mLabels.v2_gamma_type.emplace_back(0);
                }
                mLabels.v2_ancestor_pdg.emplace_back(particleTree.GetAncestorPDG(temp_trackID));
                mLabels.v2_ancestor_id.emplace_back(particleTree.GetAncestorTrackID(temp_trackID));
            }
            else
            {
                mLabels.v2_gamma_id.emplace_back(-2);
                mLabels.v2_gamma_type.emplace_back(-999);
                mLabels.v2_ancestor_pdg.emplace_back(0);
                mLabels.v2_ancestor_id.emplace_back(0);
            }
		}
        for(int j = 0; j< (int) eventArray.z2_tdc.size(); j++)
        {
            if(eventArray.z2_track_ids[j].empty() == 0 && eventArray.z2_energy[j].empty() == 0)
            {
                int temp_trackID = getLargestEnergyTrackID(eventArray.z2_track_ids[j], eventArray.z2_energy[j]);

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
                    mLabels.z2_gamma_id.emplace_back(gammas[gammaIndex].track_id);
                    mLabels.z2_gamma_type.emplace_back(gammas[gammaIndex].energy);
                }
                else
                {
                    mLabels.z2_gamma_id.emplace_back(-1);
                    mLabels.z2_gamma_type.emplace_back(0);
                }
                mLabels.z2_ancestor_pdg.emplace_back(particleTree.GetAncestorPDG(temp_trackID));
                mLabels.z2_ancestor_id.emplace_back(particleTree.GetAncestorTrackID(temp_trackID));
            }
            else
            {
                mLabels.z2_gamma_id.emplace_back(-2);
                mLabels.z2_gamma_type.emplace_back(-999);
                mLabels.z2_ancestor_pdg.emplace_back(0);
                mLabels.z2_ancestor_id.emplace_back(0);
            }
		}
        mLabelTree->Fill();
    }
}