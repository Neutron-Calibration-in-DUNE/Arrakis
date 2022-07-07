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

    void ArrayGenerator::processEvent(
        detinfo::DetectorClocksData const& clockData,
        arrakis::ParticleTree const& ParticleMaps,
        const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
        const art::ValidHandle<std::vector<raw::RawDigit>>& rawTPC
    ){
        EventArray eventArray;

        // Accquiring geometry data
        fNofAPA=fGeom->NTPC()*fGeom->Ncryostats()/2; //No. of APAs
        fChansPerAPA = fGeom->Nchannels()/fNofAPA; //No. of channels per APA

        //To get max TDC
        auto const *fDetProp = lar::providerFrom<detinfo::DetectorPropertiesService>();
        fNticks = fDetProp->NumberTimeSamples();

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

        
        unsigned int minT = 0;
        unsigned int maxT = 0;
        minT = 0;
        maxT = fNticks;
        unsigned int binT = (maxT-minT); //Bin width for TDC


        if (mcChannels.isValid() and rawTPC.isValid()){

            // Fill pointer vectors - more useful form for the raw data
            // a more usable form
            std::vector< art::Ptr<raw::RawDigit> > RawDigits;
            art::fill_ptr_vector(RawDigits, rawTPC);

            for(auto const & dptr : RawDigits) {
                const raw::RawDigit & digit = *dptr;

                // Get the channel number for this digit
                uint32_t chan = digit.Channel();
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

                //Induction Plane   
                if( fGeom->View(chan) == geo::kU){
                    for(unsigned int l=0;l<nADC_uncompPed;l++) {
                        if(uncompPed.at(l) > Threshold){

                            //Truth
                            auto const& trackIDs = truth_channel.TrackIDEs(l, l);
                            if (trackIDs.size() == 0) {
                                continue;
                            }

                            
                            if(apa == 0){
                                fEventArray.u1_tdc.emplace_back(l);
                                fEventArray.u1_channel.emplace_back(chan);
                                fEventArray.u1_adc.emplace_back( (Int_t) uncompPed.at(l));

                                std::vector<Int_t> track_ids;
                                for (size_t i = 1; i < trackIDs.size(); i++)
                                {
                                    
                                }

                                fEventArray.u1_track_ids.emplace_back(trackIDs.track_id);
                            }

                            if(apa == 1){
                                fEventArray.u1_tdc.emplace_back(l);
                                fEventArray.u1_channel.emplace_back(chan-(fNUCh+fNVCh));
                                fEventArray.u1_adc.emplace_back( (Int_t) uncompPed.at(l));
                            }

                            if(apa == 2){
                                fEventArray.u1_tdc.emplace_back(l);
                                fEventArray.u1_channel.emplace_back(chan-(fNUCh+fNVCh)*2);
                                fEventArray.u1_adc.emplace_back( (Int_t) uncompPed.at(l));
                            }

                            if(apa == 3){
                                fEventArray.u2_tdc.emplace_back(l);
                                fEventArray.u2_channel.emplace_back(chan-(fNUCh+fNVCh+fNZCh)*3);
                                fEventArray.u1_adc.emplace_back( (Int_t) uncompPed.at(l));
                            }

                            if(apa == 4){
                                fEventArray.u2_tdc.emplace_back(l);
                                fEventArray.u2_channel.emplace_back(chan-(fNUCh+fNVCh+fNZCh)*3-(fNUCh+fNVCh));
                                fEventArray.u1_adc.emplace_back( (Int_t) uncompPed.at(l));
                            }

                            if(apa == 5){
                                fEventArray.u1_tdc.emplace_back(l);
                                fEventArray.u1_channel.emplace_back(chan-(fNUCh+fNVCh+fNZCh)*3-(fNUCh+fNVCh)*2);
                                fEventArray.u1_adc.emplace_back( (Int_t) uncompPed.at(l));
                            }
                            
                            Int_t track_id = getTrackID(trackIDs, parentDaughterMap);
                            fTruthTimeChanU[apa]->Fill(truth_channel, l, particlePDGMap[track_id]);
                            //Raw data
                            fRawTimeChanU[apa]->Fill(chan, l, (Int_t) uncompPed.at(l));
                        }
                    }
                }// end of U View

                //Induction Plane   
                if( fGeom->View(chan) == geo::kV){
                    for(unsigned int l=0;l<nADC_uncompPed;l++) {
                        if(uncompPed.at(l)!=0){
                            //Truth
                            auto const& trackIDs = truth_channel.TrackIDEs(l, l);
                            if (trackIDs.size() == 0) {
                                continue;
                            }
                            Int_t track_id = getTrackID(trackIDs, parentDaughterMap);
                            fTruthTimeChanV[apa]->Fill(truth_channel, l, particlePDGMap[track_id]);
                            //Raw data
                            fRawTimeChanV[apa]->Fill(chan,l, (Int_t) uncompPed.at(l));
                        }
                    }
                }// end of V View

                //Collection Plane
                if ( fGeom->View(chan) == geo::kZ){
                    for(unsigned int l=0;l<nADC_uncompPed;l++) {
                        if(uncompPed.at(l)!=0){
                            //Truth
                            auto const& trackIDs = truth_channel.TrackIDEs(l, l);
                            if (trackIDs.size() == 0) {
                                continue;
                            }
                            Int_t track_id = getTrackID(trackIDs, parentDaughterMap);
                            fTruthTimeChanZ[apa]->Fill(truth_channel, l, particlePDGMap[track_id]);
                            //Raw data
                            fRawTimeChanZ[apa]->Fill(chan,l, (Int_t) uncompPed.at(l));
                        }
                    }
                }// end of Z View

            } //End of loop over raw digits

        } 

        //Storing the data into the root tree
        fEventArray = eventArray;
        fArrayTTree->Fill();

    }

}