/**
 * @file ArrayGenerator.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @author Junying Huang
 * @author Yashwanth Bezawada
 * @brief 
 * @version 0.1
 * @date 2022-07-07
 */
#include "ArrayGenerator.h"

namespace arrakis
{

    ArrayGenerator::ArrayGenerator(){
        fArrayTTree = fTFileService->make<TTree>("Array_Tree", "Array_Tree");
        fArrayTTree->Branch("u1_tdc", &fEventArray.u1_tdc);
        fArrayTTree->Branch("u1_channel", &fEventArray.u1_channel);
        fArrayTTree->Branch("u1_adc", &fEventArray.u1_adc);
        fArrayTTree->Branch("u1_track_ids", &fEventArray.u1_track_ids);
        fArrayTTree->Branch("u1_energy", &fEventArray.u1_energy);

        fArrayTTree->Branch("v1_tdc", &fEventArray.v1_tdc);
        fArrayTTree->Branch("v1_channel", &fEventArray.v1_channel);
        fArrayTTree->Branch("v1_adc", &fEventArray.v1_adc);
        fArrayTTree->Branch("v1_track_ids", &fEventArray.v1_track_ids);
        fArrayTTree->Branch("v1_energy", &fEventArray.v1_energy);

        fArrayTTree->Branch("z1_tdc", &fEventArray.z1_tdc);
        fArrayTTree->Branch("z1_channel", &fEventArray.z1_channel);
        fArrayTTree->Branch("z1_adc", &fEventArray.z1_adc);
        fArrayTTree->Branch("z1_track_ids", &fEventArray.z1_track_ids);
        fArrayTTree->Branch("z1_energy", &fEventArray.z1_energy);

        fArrayTTree->Branch("u2_tdc", &fEventArray.u2_tdc);
        fArrayTTree->Branch("u2_channel", &fEventArray.u2_channel);
        fArrayTTree->Branch("u2_adc", &fEventArray.u2_adc);
        fArrayTTree->Branch("u2_track_ids", &fEventArray.u2_track_ids);
        fArrayTTree->Branch("u2_energy", &fEventArray.u2_energy);

        fArrayTTree->Branch("v2_tdc", &fEventArray.v2_tdc);
        fArrayTTree->Branch("v2_channel", &fEventArray.v2_channel);
        fArrayTTree->Branch("v2_adc", &fEventArray.v2_adc);
        fArrayTTree->Branch("v2_track_ids", &fEventArray.v2_track_ids);
        fArrayTTree->Branch("v2_energy", &fEventArray.v2_energy);

        fArrayTTree->Branch("z2_tdc", &fEventArray.z2_tdc);
        fArrayTTree->Branch("z2_channel", &fEventArray.z2_channel);
        fArrayTTree->Branch("z2_adc", &fEventArray.z2_adc);
        fArrayTTree->Branch("z2_track_ids", &fEventArray.z2_track_ids);
        fArrayTTree->Branch("z2_energy", &fEventArray.z2_energy);
    }

    ArrayGenerator::~ArrayGenerator(){
    }

    void ArrayGenerator::setBoundingBoxType(std::string volumeType)
    {
        if (volumeType == "TPC" or volumeType == "tpc") { 
            fBoundingBoxType = VolumeType::TPC;
        }
        else if (volumeType == "Cryo" or volumeType == "cryo") {
            fBoundingBoxType = VolumeType::Cryostat;
        }
        else {
            fBoundingBoxType = VolumeType::World;
        }
    }

    void ArrayGenerator::ResetArrays(){
        fEventArray.u1_tdc.clear();
        fEventArray.u1_channel.clear();
        fEventArray.u1_adc.clear();
        fEventArray.u1_track_ids.clear();
        fEventArray.u1_energy.clear();
        fEventArray.v1_tdc.clear();
        fEventArray.v1_channel.clear();
        fEventArray.v1_adc.clear();
        fEventArray.v1_track_ids.clear();
        fEventArray.v1_energy.clear();
        fEventArray.z1_tdc.clear();
        fEventArray.z1_channel.clear();
        fEventArray.z1_adc.clear();
        fEventArray.z1_track_ids.clear();
        fEventArray.z1_energy.clear();
        fEventArray.u2_tdc.clear();
        fEventArray.u2_channel.clear();
        fEventArray.u2_adc.clear();
        fEventArray.u2_track_ids.clear();
        fEventArray.u2_energy.clear();
        fEventArray.v2_tdc.clear();
        fEventArray.v2_channel.clear();
        fEventArray.v2_adc.clear();
        fEventArray.v2_track_ids.clear();
        fEventArray.v2_energy.clear();
        fEventArray.z2_tdc.clear();
        fEventArray.z2_channel.clear();
        fEventArray.z2_adc.clear();
        fEventArray.z2_track_ids.clear();
        fEventArray.z2_energy.clear();
    }

    void ArrayGenerator::processEvent(
        detinfo::DetectorClocksData const& clockData,
        //arrakis::ParticleTree const& ParticleMaps,
        const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
        const art::Handle<std::vector<raw::RawDigit>>& rawTPC
    ){
        ResetArrays();

        // Accquiring geometry data
        fNofAPA=fGeom->NTPC()*fGeom->Ncryostats()/2; //No. of APAs
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
                    for(unsigned int l=0;l<nADC_uncompPed;l++) {
                        if(uncompPed.at(l) > fUPlaneThreshold){

                            auto const& trackIDs = truth_channel.TrackIDEs(l, l);
                            if (trackIDs.size() == 0) { continue; }
                                
                            std::vector<Int_t> track_ids;
                            for (size_t i = 0; i < trackIDs.size(); i++)
                            {
                                track_ids.push_back(trackIDs[i].trackID);
                            }

                            std::vector<Double_t> energy;
                            for (size_t i = 0; i < trackIDs.size(); i++)
                            {
                                energy.push_back(trackIDs[i].energy);
                            }

                            if(apa < 3){
                                fEventArray.u1_tdc.emplace_back(l);
                                fEventArray.u1_adc.emplace_back( (Int_t) uncompPed.at(l));
                                fEventArray.u1_track_ids.emplace_back(track_ids);
                                fEventArray.u1_energy.emplace_back(energy);
                            } else {
                                fEventArray.u2_tdc.emplace_back(l);
                                fEventArray.u2_adc.emplace_back( (Int_t) uncompPed.at(l));
                                fEventArray.u2_track_ids.emplace_back(track_ids);
                                fEventArray.u2_energy.emplace_back(energy);
                            }

                            if(apa == 0) {fEventArray.u1_channel.emplace_back( (Int_t) (chan) );}

                            if(apa == 1) {fEventArray.u1_channel.emplace_back( (Int_t) (chan-(fNVCh+fNZCh)) );}

                            if(apa == 2) {fEventArray.u1_channel.emplace_back( (Int_t) (chan-(fNVCh+fNZCh)*2) );}

                            if(apa == 3) {fEventArray.u2_channel.emplace_back( (Int_t) (chan-(fNUCh+fNVCh+fNZCh)*3) );}

                            if(apa == 4) {fEventArray.u2_channel.emplace_back( (Int_t) (chan-(fNUCh+fNVCh+fNZCh)*3-(fNVCh+fNZCh)) );}

                            if(apa == 5) {fEventArray.u2_channel.emplace_back( (Int_t) (chan-(fNUCh+fNVCh+fNZCh)*3-(fNVCh+fNZCh)*2) );}
                        }
                    }
                }// end of U View

                //Induction Plane V
                if( fGeom->View(chan) == geo::kV){
                    for(unsigned int l=0;l<nADC_uncompPed;l++) {
                        if(uncompPed.at(l) > fVPlaneThreshold){
                            auto const& trackIDs = truth_channel.TrackIDEs(l, l);
                            if (trackIDs.size() == 0) {
                                continue;
                            }
                            std::vector<Int_t> track_ids;
                            for (size_t i = 0; i < trackIDs.size(); i++)
                            {
                                track_ids.push_back(trackIDs[i].trackID);
                            }
                            std::vector<Double_t> energy;
                            for (size_t i = 0; i < trackIDs.size(); i++)
                            {
                                energy.push_back(trackIDs[i].energy);
                            }

                            if(apa < 3){
                                fEventArray.v1_tdc.emplace_back(l);
                                fEventArray.v1_adc.emplace_back( (Int_t) uncompPed.at(l));
                                fEventArray.v1_track_ids.emplace_back(track_ids);
                                fEventArray.v1_energy.emplace_back(energy);
                            } else {
                                fEventArray.v2_tdc.emplace_back(l);
                                fEventArray.v2_adc.emplace_back( (Int_t) uncompPed.at(l));
                                fEventArray.v2_track_ids.emplace_back(track_ids);
                                fEventArray.v2_energy.emplace_back(energy);
                            }

                            if(apa == 0) {fEventArray.v1_channel.emplace_back( (Int_t) (chan-(fNUCh)) );}

                            if(apa == 1) {fEventArray.v1_channel.emplace_back( (Int_t) (chan-(fNUCh*2+fNZCh)) );}

                            if(apa == 2) {fEventArray.v1_channel.emplace_back( (Int_t) (chan-(fNUCh*3+fNZCh*2)) );}

                            if(apa == 3) {fEventArray.v2_channel.emplace_back( (Int_t) (chan-(fNUCh+fNVCh+fNZCh)*3-(fNUCh)) );}

                            if(apa == 4) {fEventArray.v2_channel.emplace_back( (Int_t) (chan-(fNUCh+fNVCh+fNZCh)*3-(fNUCh*2+fNZCh)) );}

                            if(apa == 5) {fEventArray.v2_channel.emplace_back( (Int_t) (chan-(fNUCh+fNVCh+fNZCh)*3-(fNUCh*3+fNZCh*2)) );}
                        }
                    }
                }// end of V View

                //Collection Plane Z
                if ( fGeom->View(chan) == geo::kZ){
                    for(unsigned int l=0;l<nADC_uncompPed;l++) {
                        if(uncompPed.at(l) > fZPlaneThreshold){
                            auto const& trackIDs = truth_channel.TrackIDEs(l, l);
                            if (trackIDs.size() == 0) {
                                continue;
                            }
                            std::vector<Int_t> track_ids;
                            for (size_t i = 0; i < trackIDs.size(); i++)
                            {
                                track_ids.push_back(trackIDs[i].trackID);
                            }
                            std::vector<Double_t> energy;
                            for (size_t i = 0; i < trackIDs.size(); i++)
                            {
                                energy.push_back(trackIDs[i].energy);
                            }

                            if(apa < 3){
                                fEventArray.z1_tdc.emplace_back(l);
                                fEventArray.z1_adc.emplace_back( (Int_t) uncompPed.at(l));
                                fEventArray.z1_track_ids.emplace_back(track_ids);
                                fEventArray.z1_energy.emplace_back(energy);
                            } else {
                                fEventArray.z2_tdc.emplace_back(l);
                                fEventArray.z2_adc.emplace_back( (Int_t) uncompPed.at(l));
                                fEventArray.z2_track_ids.emplace_back(track_ids);
                                fEventArray.z2_energy.emplace_back(energy);
                            }

                            if(apa == 0) {fEventArray.z1_channel.emplace_back( (Int_t) (chan-(fNUCh+fNVCh)) );}

                            if(apa == 1) {fEventArray.z1_channel.emplace_back( (Int_t) (chan-(fNUCh+fNVCh)*2) );}

                            if(apa == 2) {fEventArray.z1_channel.emplace_back( (Int_t) (chan-(fNUCh+fNVCh)*3) );}

                            if(apa == 3) {fEventArray.z2_channel.emplace_back( (Int_t) (chan-(fNUCh+fNVCh+fNZCh)*3-(fNUCh+fNVCh)) );}

                            if(apa == 4) {fEventArray.z2_channel.emplace_back( (Int_t) (chan-(fNUCh+fNVCh+fNZCh)*3-(fNUCh+fNVCh)*2) );}

                            if(apa == 5) {fEventArray.z2_channel.emplace_back( (Int_t) (chan-(fNUCh+fNVCh+fNZCh)*3-(fNUCh+fNVCh)*3) );}
                        }
                    }
                }// end of Z View

            } //End of loop over raw digits

        } 

        //Storing the data into the root tree
        fArrayTTree->Fill();

    }

}
