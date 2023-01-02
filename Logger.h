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
#include <ctime>

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
        std::string sName;

        std::string Name() const {
            return sName;
        }

        std::string Date() const {
            return sDate;
        }

    protected:
        Logger(const std::string name);
        ~Logger() {}

    public:
        // this singleton cannot be cloned
        Logger(Logger &other) = delete;
        // singleton should also not be assignable
        void operator=(const Logger &) = delete;

        // static method that controls access to the singleton
        // instance
        static Logger* GetInstance(const std::string& name);

        std::string Preamble()
        {
            std::time_t now = std::time(0);
            std::tm *ltm = std::localtime(&now);
            std::string t = std::to_string(ltm->tm_hour) + ":" + std::to_string(ltm->tm_min) + ":" + std::to_string(ltm->tm_sec);
            std::string preamble = "[" + Date() + " " + t + "] [" + Name() + "] ";
            return preamble;
        }

        void trace(std::string status)
        {
            return Preamble() + "[trace] " + status;
        }
        
    private:
        static std::string sYear;
        static std::string sMonth;
        static std::string sDay;
        static std::string sDate;
    };    
}

#define ARRAKIS_TRACE(...)  ::arrakis::Logger::GetInstance("default")->trace(__VA_ARGS__)