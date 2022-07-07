/**
 * @file    Configuration.h
 * @author  Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief   A struct for holding LArSoft configuration parameters
 *          for the Arrakis module.
 * @version 0.1
 * @date 2022-02-15
 */
#pragma once

// art framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace arrakis
{
    struct Configuration
    {
        fhicl::Atom<bool> Generate2DArrays
        {
            fhicl::Name("Generate2DArrays"),
            fhicl::Comment("Whether to generate tdc/channel/adc arrays.")
        };
        fhicl::Atom<double> ADCThreshold
        {
            fhicl::Name("ADCThreshold"),
            fhicl::Comment("Threshold value to use for filling 2D arrays.")
        };

        fhicl::Atom<art::InputTag> LArGeantProducerLabel
        {
            fhicl::Name("LArGeantProducerLabel"),
            fhicl::Comment("Tag of the input data product for the largeant side of the simulation.")
        };

        fhicl::Atom<bool> GenerateSemanticLabels
        {
            fhicl::Name("GenerateSemanticLabels"),
            fhicl::Comment("Whether to generate semantic labels.")
        };
    };

    using Parameters = art::EDAnalyzer::Table<Configuration>;
}