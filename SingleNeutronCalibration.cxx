/**
 * @file SingleNeutronCalibration.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-08-16
 */
#include "SingleNeutronCalibration.h"

namespace arrakis
{
    SingleNeutronCalibration::SingleNeutronCalibration()
    {
        mSingleNeutronTree = mTFileService->make<TTree>("single_neutron", "single_neutron");
        mSingleNeutronTree->Branch("u1_tdc", &mSingleNeutron.u1_tdc);
        mSingleNeutronTree->Branch("u1_channel", &mSingleNeutron.u1_channel);
        mSingleNeutronTree->Branch("u1_adc", &mSingleNeutron.u1_adc);
        mSingleNeutronTree->Branch("u1_gamma_ids", &mSingleNeutron.u1_gamma_ids);
        mSingleNeutronTree->Branch("u1_gamma_energy", &mSingleNeutron.u1_gamma_energy);
        mSingleNeutronTree->Branch("u1_energy", &mSingleNeutron.u1_energy);

        mSingleNeutronTree->Branch("v1_tdc", &mSingleNeutron.v1_tdc);
        mSingleNeutronTree->Branch("v1_channel", &mSingleNeutron.v1_channel);
        mSingleNeutronTree->Branch("v1_adc", &mSingleNeutron.v1_adc);
        mSingleNeutronTree->Branch("v1_gamma_ids", &mSingleNeutron.v1_gamma_ids);
        mSingleNeutronTree->Branch("v1_gamma_energy", &mSingleNeutron.v1_gamma_energy);
        mSingleNeutronTree->Branch("v1_energy", &mSingleNeutron.v1_energy);
        
        mSingleNeutronTree->Branch("z1_tdc", &mSingleNeutron.z1_tdc);
        mSingleNeutronTree->Branch("z1_channel", &mSingleNeutron.z1_channel);
        mSingleNeutronTree->Branch("z1_adc", &mSingleNeutron.z1_adc);
        mSingleNeutronTree->Branch("z1_gamma_ids", &mSingleNeutron.z1_gamma_ids);
        mSingleNeutronTree->Branch("z1_gamma_energy", &mSingleNeutron.z1_gamma_energy);
        mSingleNeutronTree->Branch("z1_energy", &mSingleNeutron.z1_energy);

        mSingleNeutronTree->Branch("u2_tdc", &mSingleNeutron.u2_tdc);
        mSingleNeutronTree->Branch("u2_channel", &mSingleNeutron.u2_channel);
        mSingleNeutronTree->Branch("u2_adc", &mSingleNeutron.u2_adc);
        mSingleNeutronTree->Branch("u2_gamma_ids", &mSingleNeutron.u2_gamma_ids);
        mSingleNeutronTree->Branch("u2_gamma_energy", &mSingleNeutron.u2_gamma_energy);
        mSingleNeutronTree->Branch("u2_energy", &mSingleNeutron.u2_energy);

        mSingleNeutronTree->Branch("v2_tdc", &mSingleNeutron.v2_tdc);
        mSingleNeutronTree->Branch("v2_channel", &mSingleNeutron.v2_channel);
        mSingleNeutronTree->Branch("v2_adc", &mSingleNeutron.v2_adc);
        mSingleNeutronTree->Branch("v2_gamma_ids", &mSingleNeutron.v2_gamma_ids);
        mSingleNeutronTree->Branch("v2_gamma_energy", &mSingleNeutron.v2_gamma_energy);
        mSingleNeutronTree->Branch("v2_energy", &mSingleNeutron.v2_energy);

        mSingleNeutronTree->Branch("z2_tdc", &mSingleNeutron.z2_tdc);
        mSingleNeutronTree->Branch("z2_channel", &mSingleNeutron.z2_channel);
        mSingleNeutronTree->Branch("z2_adc", &mSingleNeutron.z2_adc);
        mSingleNeutronTree->Branch("z2_gamma_ids", &mSingleNeutron.z2_gamma_ids);
        mSingleNeutronTree->Branch("z2_gamma_energy", &mSingleNeutron.z2_gamma_energy);
        mSingleNeutronTree->Branch("z2_energy", &mSingleNeutron.z2_energy);
    }

    SingleNeutronCalibration::~SingleNeutronCalibration()
    {
    }

    void SingleNeutronCalibration::ResetSingleNeutron()
    {
        mSingleNeutron.u1_tdc.clear();
        mSingleNeutron.u1_channel.clear();
        mSingleNeutron.u1_adc.clear();
        mSingleNeutron.u1_gamma_ids.clear();
        mSingleNeutron.u1_gamma_energy.clear();
        mSingleNeutron.u1_energy.clear();

        mSingleNeutron.v1_tdc.clear();
        mSingleNeutron.v1_channel.clear();
        mSingleNeutron.v1_adc.clear();
        mSingleNeutron.v1_gamma_ids.clear();
        mSingleNeutron.v1_gamma_energy.clear();
        mSingleNeutron.v1_energy.clear();
        
        mSingleNeutron.z1_tdc.clear();
        mSingleNeutron.z1_channel.clear();
        mSingleNeutron.z1_adc.clear();
        mSingleNeutron.z1_gamma_ids.clear();
        mSingleNeutron.z1_gamma_energy.clear();
        mSingleNeutron.z1_energy.clear();
        
        mSingleNeutron.u2_tdc.clear();
        mSingleNeutron.u2_channel.clear();
        mSingleNeutron.u2_adc.clear();
        mSingleNeutron.u2_gamma_ids.clear();
        mSingleNeutron.u2_gamma_energy.clear();
        mSingleNeutron.u2_energy.clear();
        
        mSingleNeutron.v2_tdc.clear();
        mSingleNeutron.v2_channel.clear();
        mSingleNeutron.v2_adc.clear();
        mSingleNeutron.v2_gamma_ids.clear();
        mSingleNeutron.v2_gamma_energy.clear();
        mSingleNeutron.v2_energy.clear();
        
        mSingleNeutron.z2_tdc.clear();
        mSingleNeutron.z2_channel.clear();
        mSingleNeutron.z2_adc.clear();
        mSingleNeutron.z2_gamma_ids.clear();
        mSingleNeutron.z2_gamma_energy.clear();
        mSingleNeutron.z2_energy.clear();
    }

    void SingleNeutronCalibration::processEvent(
        ParticleTree particleTree,
        // arrakis::ParticleTree const& ParticleMaps,
        // const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
        const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
        const art::Handle<std::vector<raw::RawDigit>>& rawTPC
        // const art::ValidHandle<std::vector<sim::SimEnergyDeposit>>& mcEnergyDeposits
    )
    {
        ResetSingleNeutron();

        // Accquiring geometry data
        fNofAPA = 6;//fGeom->NTPC()*fGeom->Ncryostats()/2; //No. of APAs
        fChansPerAPA = fGeom->Nchannels()/fNofAPA; //No. of channels per APA

        // taken from dune35t module a way to organise the channel mapping:
        // loop through channels in the first APA to find the channel boundaries for each view
        // will adjust for desired APA after
        fUChanMin = 0;
        fZChanMax = fChansPerAPA - 1;
        for ( unsigned int c = fUChanMin + 1; c < fZChanMax; c++ ){
            if ( fGeom->View(c) == geo::kV && fGeom->View(c-1) == geo::kU ){
                fVChanMin = c;
                fUChanMax = c - 1;
            }
            if ( fGeom->View(c) == geo::kZ && fGeom->View(c-1) == geo::kV ){
                fZChanMin = c;
                fVChanMax = c-1;
            }
        }

        //Number of channels in each view
        fNUCh=fUChanMax-fUChanMin+1; //U
        fNVCh=fVChanMax-fVChanMin+1; //V
        fNZCh=fZChanMax-fZChanMin+1; //Z (collection plane)

        if (mcChannels.isValid() and rawTPC.isValid()){

            // Fill pointer vectors - more useful form for the raw data
            // a more usable form
            std::vector< art::Ptr<raw::RawDigit> > RawDigits;
            art::fill_ptr_vector(RawDigits, rawTPC);

            for(auto const & dptr : RawDigits) {
                const raw::RawDigit & digit = *dptr;

                // Get the channel number for this digit
                auto chan = digit.Channel();
                // number of samples in uncompressed ADC
                int nSamples = digit.Samples();
                unsigned int apa = std::floor( chan/fChansPerAPA );	  
                int pedestal = (int)digit.GetPedestal();
                
                std::vector<short> uncompressed(nSamples);
                // with pedestal	  
                raw::Uncompress(digit.ADCs(), uncompressed, pedestal, digit.Compression());
                // subtract pedestals
                std::vector<short> uncompPed(nSamples);
                for (int i=0; i<nSamples; i++) uncompPed.at(i)=uncompressed.at(i)-pedestal;
                
                // number of ADC uncompressed without pedestal
                nADC_uncompPed=uncompPed.size();	  

                //Truth Channel
                auto truth_channel = mcChannels->at(chan); 

                //Induction Plane U
                if( fGeom->View(chan) == geo::kU){
                    //Loop over TDC
                    for(unsigned int l=0;l<nADC_uncompPed;l++) {
                        //Checking for the ADC threshold
                        if(std::abs(uncompPed.at(l)) > fUPlaneThreshold){

                            auto const& trackIDs = truth_channel.TrackIDEs(l, l);
                            if (trackIDs.size() == 0) { continue; }

                            for (size_t i = 0; i < trackIDs.size(); i++){
                                Double_t eFrac = trackIDs[i].energyFrac;
                                Double_t energy = trackIDs[i].energy;
                                Int_t ancestor_pdg = particleTree.GetAncestorPDG( trackIDs[i].trackID );
                                Int_t gamma_id = 0;
                                Double_t gamma_energy = 0;
                                if (ancestor_pdg == 22)
                                {
                                    gamma_id = particleTree.GetAncestorTrackID( trackIDs[i].trackID );
                                    gamma_energy = particleTree.GetAncestorEnergy( trackIDs[i].trackID );
                                } else {
                                    gamma_id = -1;
                                    gamma_energy = -1;
                                }

                                int a = std::floor((apa+1)/2);
                                
                                if(apa == 0 || apa == 2 || apa == 4){
                                    mSingleNeutron.u1_tdc.emplace_back(l);
                                    mSingleNeutron.u1_channel.emplace_back( (Int_t) (chan-(fNUCh*a+(fNVCh+fNZCh)*apa) ));
                                    mSingleNeutron.u1_adc.emplace_back( (Int_t) (std::abs(uncompPed.at(l))*eFrac) );
                                    mSingleNeutron.u1_gamma_ids.emplace_back(gamma_id);
                                    mSingleNeutron.u1_gamma_energy.emplace_back(gamma_energy);
                                    mSingleNeutron.u1_energy.emplace_back(energy);
                                } else {
                                    mSingleNeutron.u2_tdc.emplace_back(l);
                                    mSingleNeutron.u2_channel.emplace_back( (Int_t) (chan-(fNUCh*a+(fNVCh+fNZCh)*apa) ));
                                    mSingleNeutron.u2_adc.emplace_back( (Int_t) (std::abs(uncompPed.at(l))*eFrac) );
                                    mSingleNeutron.u2_gamma_ids.emplace_back(gamma_id);
                                    mSingleNeutron.u2_gamma_energy.emplace_back(gamma_energy);
                                    mSingleNeutron.u2_energy.emplace_back(energy);
                                }

                                // if(apa == 0) {mSingleNeutron.u1_channel.emplace_back( (Int_t) (chan) );}

                                // if(apa == 1) {mSingleNeutron.u2_channel.emplace_back( (Int_t) (chan-(fNUCh+fNVCh+fNZCh)) );}

                                // if(apa == 2) {mSingleNeutron.u1_channel.emplace_back( (Int_t) (chan-(fNUCh+(fNVCh+fNZCh)*2)) );}

                                // if(apa == 3) {mSingleNeutron.u2_channel.emplace_back( (Int_t) (chan-(fNUCh*2+(fNVCh+fNZCh)*3)) );}

                                // if(apa == 4) {mSingleNeutron.u1_channel.emplace_back( (Int_t) (chan-(fNUCh*2+(fNVCh+fNZCh)*4)) );}

                                // if(apa == 5) {mSingleNeutron.u2_channel.emplace_back( (Int_t) (chan-(fNUCh*3+(fNVCh+fNZCh)*5) );}

                            }
                        }
                    }
                }// end of U View

                //Induction Plane V
                if( fGeom->View(chan) == geo::kV){
                    for(unsigned int l=0;l<nADC_uncompPed;l++) {
                        if(std::abs(uncompPed.at(l)) > fVPlaneThreshold){
                            auto const& trackIDs = truth_channel.TrackIDEs(l, l);
                            if (trackIDs.size() == 0) {
                                continue;
                            }

                            for (size_t i = 0; i < trackIDs.size(); i++){
                                Double_t eFrac = trackIDs[i].energyFrac;
                                Double_t energy = trackIDs[i].energy;
                                Int_t ancestor_pdg = particleTree.GetAncestorPDG( trackIDs[i].trackID );
                                Int_t gamma_id = 0;
                                Double_t gamma_energy = 0;
                                if (ancestor_pdg == 22)
                                {
                                    gamma_id = particleTree.GetAncestorTrackID( trackIDs[i].trackID );
                                    gamma_energy = particleTree.GetAncestorEnergy( trackIDs[i].trackID );
                                } else {
                                    gamma_id = -1;
                                    gamma_energy = -1;
                                }

                                int a = std::floor((apa+1)/2);
                    
                                if(apa == 0 || apa == 2 || apa == 4){
                                    mSingleNeutron.v1_tdc.emplace_back(l);
                                    mSingleNeutron.v1_channel.emplace_back( (Int_t) (chan-(fNUCh*(apa+1)+fNVCh*a+fNZCh*apa) ));
                                    mSingleNeutron.v1_adc.emplace_back( (Int_t) (std::abs(uncompPed.at(l))*eFrac) );
                                    mSingleNeutron.v1_gamma_ids.emplace_back(gamma_id);
                                    mSingleNeutron.v1_gamma_energy.emplace_back(gamma_energy);
                                    mSingleNeutron.v1_energy.emplace_back(energy);
                                } else {
                                    mSingleNeutron.v2_tdc.emplace_back(l);
                                    mSingleNeutron.v2_channel.emplace_back( (Int_t) (chan-(fNUCh*(apa+1)+fNVCh*a+fNZCh*apa) ));
                                    mSingleNeutron.v2_adc.emplace_back( (Int_t) (std::abs(uncompPed.at(l))*eFrac) );
                                    mSingleNeutron.v2_gamma_ids.emplace_back(gamma_id);
                                    mSingleNeutron.v2_gamma_energy.emplace_back(gamma_energy);
                                    mSingleNeutron.v2_energy.emplace_back(energy);
                                }

                                // if(apa == 0) {mSingleNeutron.v1_channel.emplace_back( (Int_t) (chan-(fNUCh)) );}

                                // if(apa == 1) {mSingleNeutron.v2_channel.emplace_back( (Int_t) (chan-(fNUCh*2+fNVCh+fNZCh)) );}

                                // if(apa == 2) {mSingleNeutron.v1_channel.emplace_back( (Int_t) (chan-(fNUCh*3+fNVCh+fNZCh*2)) );}

                                // if(apa == 3) {mSingleNeutron.v2_channel.emplace_back( (Int_t) (chan-(fNUCh*4+fNVCh*2+fNZCh*3)) );}

                                // if(apa == 4) {mSingleNeutron.v1_channel.emplace_back( (Int_t) (chan-(fNUCh*5+fNVCh*2+fNZCh*4)) );}

                                // if(apa == 5) {mSingleNeutron.v2_channel.emplace_back( (Int_t) (chan-(fNUCh*6+fNVCh*3+fNZCh*5)) );}

                            }
                        }
                    }
                }// end of V View

                //Collection Plane Z
                if ( fGeom->View(chan) == geo::kZ){
                    for(unsigned int l=0;l<nADC_uncompPed;l++) {
                        if(std::abs(uncompPed.at(l)) > fZPlaneThreshold){
                            auto const& trackIDs = truth_channel.TrackIDEs(l, l);
                            if (trackIDs.size() == 0) {
                                continue;
                            }

                            for (size_t i = 0; i < trackIDs.size(); i++){
                                Double_t eFrac = trackIDs[i].energyFrac;
                                Double_t energy = trackIDs[i].energy;
                                Int_t ancestor_pdg = particleTree.GetAncestorPDG( trackIDs[i].trackID );
                                Int_t gamma_id = 0;
                                Double_t gamma_energy = 0;
                                if (ancestor_pdg == 22)
                                {
                                    gamma_id = particleTree.GetAncestorTrackID( trackIDs[i].trackID );
                                    gamma_energy = particleTree.GetAncestorEnergy( trackIDs[i].trackID );
                                } else {
                                    gamma_id = -1;
                                    gamma_energy = -1;
                                }

                                int a = std::floor((apa+1)/2);

                                if(apa == 0 || apa == 2 || apa == 4){
                                    mSingleNeutron.z1_tdc.emplace_back(l);
                                    mSingleNeutron.z1_channel.emplace_back( (Int_t) (chan-((fNUCh+fNVCh)*(apa+1)+fNZCh*a) ));
                                    mSingleNeutron.z1_adc.emplace_back( (Int_t) (std::abs(uncompPed.at(l))*eFrac) );
                                    mSingleNeutron.z1_gamma_ids.emplace_back(gamma_id);
                                    mSingleNeutron.z1_gamma_energy.emplace_back(gamma_energy);
                                    mSingleNeutron.z1_energy.emplace_back(energy);
                                } else {
                                    mSingleNeutron.z2_tdc.emplace_back(l);
                                    mSingleNeutron.z2_channel.emplace_back( (Int_t) (chan-((fNUCh+fNVCh)*(apa+1)+fNZCh*a) ));
                                    mSingleNeutron.z2_adc.emplace_back( (Int_t) (std::abs(uncompPed.at(l))*eFrac) );
                                    mSingleNeutron.z2_gamma_ids.emplace_back(gamma_id);
                                    mSingleNeutron.z2_gamma_energy.emplace_back(gamma_energy);
                                    mSingleNeutron.z2_energy.emplace_back(energy);
                                }

                                // if(apa == 0) {mSingleNeutron.z1_channel.emplace_back( (Int_t) (chan-(fNUCh+fNVCh)) );}

                                // if(apa == 1) {mSingleNeutron.z2_channel.emplace_back( (Int_t) (chan-((fNUCh+fNVCh)*2+fNZCh)) );}

                                // if(apa == 2) {mSingleNeutron.z1_channel.emplace_back( (Int_t) (chan-((fNUCh+fNVCh)*3+fNZCh)) );}

                                // if(apa == 3) {mSingleNeutron.z2_channel.emplace_back( (Int_t) (chan-((fNUCh+fNVCh)*4+fNZCh*2)) );}

                                // if(apa == 4) {mSingleNeutron.z1_channel.emplace_back( (Int_t) (chan-((fNUCh+fNVCh)*5+fNZCh*2)) );}

                                // if(apa == 5) {mSingleNeutron.z2_channel.emplace_back( (Int_t) (chan-((fNUCh+fNVCh)*6+fNZCh*3)) );}

                            }
                        }
                    }
                }// end of Z View
            } //End of loop over raw digits
        }
        mSingleNeutronTree->Fill();
    }
}

// mSingleNeutron.edep_x.clear();
// mSingleNeutron.edep_y.clear();
// mSingleNeutron.edep_z.clear();
// mSingleNeutron.edep_energy.clear();
// mSingleNeutron.edep_num_electrons.clear();
// mSingleNeutron.edep_gamma_ids.clear();
// mSingleNeutron.edep_total_energy = 0.0;
// mSingleNeutron.complete_apa = true;
// mSingleNeutron.complete_capture = true;



// mSingleNeutronTree->Branch("edep_x", &mSingleNeutron.edep_x);
// mSingleNeutronTree->Branch("edep_y", &mSingleNeutron.edep_y);
// mSingleNeutronTree->Branch("edep_z", &mSingleNeutron.edep_z);
// mSingleNeutronTree->Branch("edep_energy", &mSingleNeutron.edep_energy);
// mSingleNeutronTree->Branch("edep_num_electrons", &mSingleNeutron.edep_num_electrons);
// mSingleNeutronTree->Branch("edep_gamma_ids", &mSingleNeutron.edep_gamma_ids);
// mSingleNeutronTree->Branch("edep_total_energy", &mSingleNeutron.edep_total_energy);
// mSingleNeutronTree->Branch("complete_apa", &mSingleNeutron.complete_apa);
// mSingleNeutronTree->Branch("complete_capture", &mSingleNeutron.complete_capture);



// bool complete_apa = true;
// bool positive = false;
// bool negative = false;
// for (auto energyDeposit : *mcEnergyDeposits)
// {
//     if(
//         particleTree.GetAncestorPDG(energyDeposit.TrackID()) == 2112 &&
//         particleTree.GetParentPDG(energyDeposit.TrackID()) != 2112 && 
//         particleTree.GetPDGCode(energyDeposit.TrackID()) != 2112 &&
//         particleTree.GetPDGCode(energyDeposit.TrackID()) == 11
//     )
//     {
//         mSingleNeutron.edep_x.emplace_back(energyDeposit.StartX());
//         if (energyDeposit.StartX() <= 0.0) {
//             negative = true;
//         }
//         else {
//             positive = true;
//         }
//         mSingleNeutron.edep_y.emplace_back(energyDeposit.StartY());
//         mSingleNeutron.edep_z.emplace_back(energyDeposit.StartZ());
//         mSingleNeutron.edep_energy.emplace_back(energyDeposit.Energy());
//         mSingleNeutron.edep_num_electrons.emplace_back(energyDeposit.NumElectrons());
//         mSingleNeutron.edep_total_energy += energyDeposit.Energy();

//         Int_t pdg = particleTree.GetPDGCode(energyDeposit.TrackID());
//         Int_t track_energy = energyDeposit.TrackID();
//         while (pdg != 22)
//         {
//             if (pdg == 2112) {
//                 break;
//             }                    
//             Int_t temp_track_energy = particleTree.GetParentTrackID(track_energy);
//             pdg = particleTree.GetPDGCode(temp_track_energy);
//             track_energy = temp_track_energy;
//         }
//         mSingleNeutron.edep_gamma_ids.emplace_back(track_energy);
//     }
// }
// if (positive && negative) {
//     complete_apa = false;
// }
// mSingleNeutron.complete_apa = complete_apa;
// if (round(mSingleNeutron.edep_total_energy) != 6.1) { 
//     mSingleNeutron.complete_capture = false;
// }