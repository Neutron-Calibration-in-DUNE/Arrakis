/**
 * @file Logger.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-07-07
 */
#include "Logger.h"

namespace arrakis 
{
    Logger* Logger::sInstance{nullptr};
    std::mutex Logger::sMutex;

    Logger *Logger::GetInstance()
    {
        std::lock_guard<std::mutex> lock(sMutex);
        if (sInstance == nullptr)
        {
            sInstance = new Logger();
        }
        return sInstance;
    }

    Logger::Logger()
    {
    }

    
}