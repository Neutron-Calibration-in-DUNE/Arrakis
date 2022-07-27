/**
 * @file EvtLvlNeutronInfo.cxx
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-07-26
 */
#include "EvtLvlNeutronInfo.h"

namespace arrakis
{
    EvtLvlNeutronInfo::EvtLvlNeutronInfo()
    {
        mEventStatisticsTree = mTFileService->make<TTree>("event_statistics", "event_statistics");

        mEventStatisticsTree->Branch("total_neutrons_captured", &mEventStatistics.total_neutrons_captured);

        mEventStatisticsTree->Branch("u_summed_adc", &mEventStatistics.u_summed_adc);
        mEventStatisticsTree->Branch("u_total_summed_adc", &mEventStatistics.u_total_summed_adc);
        mEventStatisticsTree->Branch("u_total_summed_energy", &mEventStatistics.u_total_summed_energy);

        mEventStatisticsTree->Branch("v_summed_adc", &mEventStatistics.v_summed_adc);
        mEventStatisticsTree->Branch("v_total_summed_adc", &mEventStatistics.v_total_summed_adc);
        mEventStatisticsTree->Branch("v_total_summed_energy", &mEventStatistics.v_total_summed_energy);

        mEventStatisticsTree->Branch("z_summed_adc", &mEventStatistics.z_summed_adc);
        mEventStatisticsTree->Branch("z_total_summed_adc", &mEventStatistics.z_total_summed_adc);
        mEventStatisticsTree->Branch("z_total_summed_energy", &mEventStatistics.z_total_summed_energy);
    }

    EvtLvlNeutronInfo::~EvtLvlNeutronInfo()
    {}

    void EvtLvlNeutronInfo::processEvent(
        detinfo::DetectorClocksData const& clockData,
        const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
        const art::Handle<std::vector<raw::RawDigit>>& rawTPC
    )
    {
        if (mcChannels.isValid() and mcParticles.isValid() and rawTPC.isValid())
        {
            // Fill pointer vectors - more useful form for the raw data
            // a more usable form
            std::vector< art::Ptr<raw::RawDigit> > RawDigits;
            art::fill_ptr_vector(RawDigits, rawTPC);
            
            // Struct to store information in
            EventStatistics event_statistics;

            Int_t neutrons_captured = 0;

            std::vector<Int_t> u_adc(fClockTicks, 0);
            Int_t u_summed_adc = 0;
            Double_t u_summed_energy = 0;

            std::vector<Int_t> v_adc(fClockTicks, 0);
            Int_t v_summed_adc = 0;
            Double_t v_summed_energy = 0;

            std::vector<Int_t> z_adc(fClockTicks, 0);
            Int_t z_summed_adc = 0;
            Double_t z_summed_energy = 0;

            /////////// Calculating number of neutron captures ////////////////////////////
            for (auto particle : *mcParticles) {
                if (particle.PdgCode() == 2112) //If the particle is neutron
                {
                    if (particle.EndProcess() == "nCapture") //If the neutron is captured
                    {
                        DetectorVolume ending_volume = fGeometry->getVolume(
                            particle.EndX(), particle.EndY(), particle.EndZ()
                        );
                        if (ending_volume.material_name == "LAr" and ending_volume.volume_type == 2) { //If the neutron is captured in the TPC
                            neutrons_captured += 1;
                        }
                    }
                }
            }
            
            /////////// Calculating summed ADC and summed energy //////////////////////////
            for(auto const & dptr : RawDigits) {
                const raw::RawDigit & digit = *dptr;

                // Get the channel number for this digit
                auto chan = digit.Channel();
                // number of samples in uncompressed ADC
                int nSamples = digit.Samples();
                // Pedestal level (ADC counts) 
                int pedestal = (int)digit.GetPedestal();
                // To store uncompressed ADC samples
                std::vector<short> uncompressed(nSamples);
                // uncompress the raw data buffer
                raw::Uncompress(digit.ADCs(), uncompressed, pedestal, digit.Compression());
                // subtract pedestals from the uncompressed ADCs
                std::vector<short> uncompPed(nSamples);
                for (int i=0; i<nSamples; i++) uncompPed.at(i)=uncompressed.at(i)-pedestal;
                //Truth Channel
                auto truth_channel = mcChannels->at(chan);

                // number of ADC uncompressed without pedestal
                clock_ticks=uncompPed.size();
                
                /////////////////////////////////////////// U Channel ///////////////////////////////////////
                if( fGeom->View(chan) == geo::kU) {
                    //loop over time ticks
                    for(unsigned int l=0;l<clock_ticks;l++) {
                        if(uncompPed.at(l) > fUPlaneThreshold){
                            auto const& trackIDs = truth_channel.TrackIDEs(l, l);
                            if (trackIDs.size() == 0) { continue; }

                            for (size_t i = 0; i < trackIDs.size(); i++)
                            {
                                u_summed_energy += trackIDs[i].energy;
                            }

                            u_adc.at( (int) l) += uncompPed.at(l);
                            u_summed_adc += (Int_t) uncompPed.at(l);
                        }
                    }
                }// end of U View

                /////////////////////////////////////////// V Channel ///////////////////////////////////////
                if( fGeom->View(chan) == geo::kV) {
                    //loop over time ticks
                    for(unsigned int l=0;l<clock_ticks;l++) {
                        if(uncompPed.at(l) > fVPlaneThreshold){
                            auto const& trackIDs = truth_channel.TrackIDEs(l, l);
                            if (trackIDs.size() == 0) { continue; }

                            for (size_t i = 0; i < trackIDs.size(); i++)
                            {
                                v_summed_energy += trackIDs[i].energy;
                            }
                            
                            v_adc.at( (int) l) += uncompPed.at(l);
                            v_summed_adc += (Int_t) uncompPed.at(l);
                        }
                    }
                }// end of V View

                /////////////////////////////////////////// Z Channel ///////////////////////////////////////
                if( fGeom->View(chan) == geo::kZ) {
                    //loop over time ticks
                    for(unsigned int l=0;l<clock_ticks;l++) {
                        if(uncompPed.at(l) > fZPlaneThreshold){
                            auto const& trackIDs = truth_channel.TrackIDEs(l, l);
                            if (trackIDs.size() == 0) { continue; }

                            for (size_t i = 0; i < trackIDs.size(); i++)
                            {
                                z_summed_energy += trackIDs[i].energy;
                            }
                            
                            z_adc.at( (int) l) += uncompPed.at(l);
                            z_summed_adc += (Int_t) uncompPed.at(l);
                        }
                    }
                }// end of Z View
            }

            //////// Filling the tree
            event_statistics.total_neutrons_captured = neutrons_captured;

            event_statistics.u_summed_adc = u_adc;
            event_statistics.u_total_summed_adc = u_summed_adc;
            event_statistics.u_total_summed_energy = u_summed_energy;

            event_statistics.v_summed_adc = v_adc;
            event_statistics.v_total_summed_adc = v_summed_adc;
            event_statistics.v_total_summed_energy = v_summed_energy;

            event_statistics.z_summed_adc = z_adc;
            event_statistics.z_total_summed_adc = z_summed_adc;
            event_statistics.z_total_summed_energy = z_summed_energy;

            mEventStatistics = event_statistics;
            mEventStatisticsTree->Fill();
        }
    }
}