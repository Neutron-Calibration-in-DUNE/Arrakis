/**
 * @file Generators.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-12/21
 */
#include "Generators.h"

namespace arrakis
{
    Generators::Generators()
    {
    }
    Generators::~Generators()
    {
    }

    void Generators::ProcessEvent(
        std::vector<art::InputTag> inputTags,
        std::vector<std::string> labels,
        std::vector<art::ValidHandle<simb::MCTruth>> mcTruth
    )
    {
        for(size_t ii = 0; ii < labels.size(); ii++)
        {
            mGenerators.emplace_back(
                Generator(inputTags[ii], labels[ii], mcTruth[ii])
            );
        }
    }
}