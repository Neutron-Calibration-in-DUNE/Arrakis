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
            FG_DEFAULT  = 39,
            FG_BLACK    = 30,
            FG_RED      = 31,
            FG_GREEN    = 32,
            FG_YELLOW   = 33,
            FG_BLUE     = 34,
            FG_MAGENTA  = 35,
            FG_CYAN     = 36,
            FG_LIGHT_GRAY   = 37,
            FG_DARK_GRAY    = 90,
            FG_LIGHT_RED    = 91,
            FG_LIGHT_GREEN  = 92,
            FG_LIGHT_YELLOW = 93,
            FG_LIGHT_BLUE   = 94,
            FG_LIGHT_MAGENTA = 95, 
            FG_LIGHT_CYAN   = 96, 
            FG_WHITE        = 97,
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
    static Color::Modifier def = Color::Modifier(Color::FG_DEFAULT);
    static Color::Modifier black = Color::Modifier(Color::FG_BLACK);
    static Color::Modifier red = Color::Modifier(Color::FG_RED);
    static Color::Modifier green = Color::Modifier(Color::FG_GREEN);
    static Color::Modifier yellow = Color::Modifier(Color::FG_YELLOW);
    static Color::Modifier blue = Color::Modifier(Color::FG_BLUE);
    static Color::Modifier magenta = Color::Modifier(Color::FG_MAGENTA);
    static Color::Modifier cyan = Color::Modifier(Color::FG_CYAN);
    static Color::Modifier light_gray = Color::Modifier(Color::FG_LIGHT_GRAY);
    static Color::Modifier dark_gray = Color::Modifier(Color::FG_DARK_GRAY);
    

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
            return "[" + Date() + " " + t + "]";
        }

        void WriteToFile(std::string status)
        {
            std::ofstream file;
            file.open(sOutputFileName, std::ios::out | std::ios::app);
            if(file.fail()) {
                throw std::ios_base::failure(std::strerror(errno));
            }
            file.exceptions(file.exceptions() | std::ios::failbit | std::ifstream::badbit);
            file << status << std::endl;
        }

        void trace(std::string status)
        {
            std::cout << Preamble() + " [" << green << Name() << def << "] ";
            std::cout << "[" << light_gray << "trace" << def << "] " + status << std::endl;
            std::string file_output = Preamble() + " [" + Name() + "] [trace] " + status;
            WriteToFile(file_output);
        }
        
    private:
        std::string sYear;
        std::string sMonth;
        std::string sDay;
        std::string sDate;
        std::string sOutputFileName;
    };    
}

#define ARRAKIS_TRACE(...)  ::arrakis::Logger::GetInstance("default")->trace(__VA_ARGS__)