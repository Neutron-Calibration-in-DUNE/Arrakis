/**
 * @file    Configuration.h
 * @author  Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief   A struct for holding LArSoft configuration parameters
 *          for the Arrakis module.
 * @version 0.1
 * @date 2022-02-15
 */
#pragma once
#include <string>
#include <vector>
#include <memory>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <mutex>

namespace arrakis
{
    /**
     * @brief Singleton class for logging information.
     * 
     */
    class Logger
    {
    private:
        static Logger * sInstance;
        static std::mutex sMutex;

    protected:
        Logger();
        ~Logger() {}

    public:
        // this singleton cannot be cloned
        Logger(Logger &other) = delete;
        // singleton should also not be assignable
        void operator=(const Logger &) = delete;

        // static method that controls access to the singleton
        // instance
        static Logger* GetInstance();

        static void trace(std::string status)
        {
            std::cout << status << std::endl;
        }
        
    private:
        
    };    
}

#define ARRAKIS_TRACE(...)  ::arrakis::Logger::GetInstance()->trace(__VA_ARGS__)