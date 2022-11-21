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

        mSingleNeutronTree->Branch("tdc", &mSingleNeutron.tdc);
        mSingleNeutronTree->Branch("channel", &mSingleNeutron.channel);
        mSingleNeutronTree->Branch("adc", &mSingleNeutron.adc);
        // mSingleNeutronTree->Branch("gamma_ids", &mSingleNeutron.gamma_ids);
        // mSingleNeutronTree->Branch("gamma_energy", &mSingleNeutron.gamma_energy);
        mSingleNeutronTree->Branch("energy", &mSingleNeutron.energy);
    }

    SingleNeutronCalibration::~SingleNeutronCalibration()
    {
    }

    void SingleNeutronCalibration::ResetSingleNeutron()
    {
        mSingleNeutron.tdc.clear();
        mSingleNeutron.channel.clear();
        mSingleNeutron.adc.clear();
        // mSingleNeutron.gamma_ids.clear();
        // mSingleNeutron.gamma_energy.clear();
        mSingleNeutron.energy.clear();
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
                // unsigned int apa = std::floor( chan/fChansPerAPA );	  
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

                //Setting threshold
                Double_t threshold = 0;
                if( fGeom->View(chan) == geo::kU){
                    threshold = fUPlaneThreshold;
                }
                if( fGeom->View(chan) == geo::kV){
                    threshold = fVPlaneThreshold;
                }
                if( fGeom->View(chan) == geo::kZ){
                    threshold = fZPlaneThreshold;
                }
                
                for(unsigned int l=0;l<nADC_uncompPed;l++) {
                    if(std::abs(uncompPed.at(l)) > threshold){

                        auto const& trackIDs = truth_channel.TrackIDEs(l, l);
                        if (trackIDs.size() == 0) { continue; }

                        for (size_t i = 0; i < trackIDs.size(); i++){
                            Double_t eFrac = trackIDs[i].energyFrac;
                            Double_t energy = trackIDs[i].energy;
                            // Int_t ancestor_pdg = particleTree.GetAncestorPDG( trackIDs[i].trackID );
                            // Int_t gamma_id = 0;
                            // Double_t gamma_energy = 0;
                            // if (ancestor_pdg == 22)
                            // {
                            //     gamma_id = particleTree.GetAncestorTrackID( trackIDs[i].trackID );
                            //     gamma_energy = particleTree.GetAncestorEnergy( trackIDs[i].trackID );
                            // } else {
                            //     gamma_id = -1;
                            //     gamma_energy = -1;
                            // }

                            mSingleNeutron.tdc.emplace_back(l);
                            mSingleNeutron.channel.emplace_back( (Int_t) (chan));
                            mSingleNeutron.adc.emplace_back( (Int_t) (std::abs(uncompPed.at(l))*eFrac) );
                            // mSingleNeutron.gamma_ids.emplace_back(gamma_id);
                            // mSingleNeutron.gamma_energy.emplace_back(gamma_energy);
                            mSingleNeutron.energy.emplace_back(energy);
                        }
                    }
                }   
            } //End of loop over raw digits
        }
        mSingleNeutronTree->Fill();
    }
}