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
    namespace Color 
    {
        enum Code {
            FG_RED      = 31,
            FG_GREEN    = 32,
            FG_BLUE     = 34,
            FG_DEFAULT  = 39,
            BG_RED      = 41,
            BG_GREEN    = 42,
            BG_BLUE     = 44,
            BG_DEFAULT  = 49
        };
        class Modifier {
            Code code;
        public:
            Modifier(Code pCode) : code(pCode) {}
            friend std::ostream&
            operator<<(std::ostream& os, const Modifier& mod) {
                return os << "\033[" << mod.code << "m";
            }
        };
    }
    static Color::Modifier red = Color::Modifier(Color::FG_RED);
    static Color::Modifier def = Color::Modifier(Color::FG_DEFAULT);
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
            std::string t = std::to_string(ltm->tm_hour);
            t += ":" + std::to_string(ltm->tm_min);
            t += ":" + std::to_string(ltm->tm_sec);
            return "[" + Date() + " " + t + "] [" + Name() + "] ";
        }

        void trace(std::string status)
        {
            std::cout << Preamble() + "[" << red << "trace" << def << "] " + status << std::endl;
        }
        
    private:
        std::string sYear;
        std::string sMonth;
        std::string sDay;
        std::string sDate;
    };    
}

#define ARRAKIS_TRACE(...)  ::arrakis::Logger::GetInstance("default")->trace(__VA_ARGS__)