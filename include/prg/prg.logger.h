#ifndef PRG_LOGGER
#define PRG_LOGGER

/**
* @file Отображение памяти в файл и чтение памяти из файла.
*/

#include <iostream>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <fstream>
#include <string>
#include "prg.string.h"

namespace prg
{
    struct Logger : public std::string
    {
        std::string filename;
        bool print_flag;

        Logger(bool flag, const std::string &filename="") : filename(filename), print_flag(flag) {}
        ~Logger() { if (print_flag) flush(); }

        void insert(const char *format) { this->append(prg::to_string(format)); }

        template<typename T, typename... Args>
        void insert(const char *format, T value, Args... args)
            { this->append(prg::to_string(format, value, args...)); }

        void flush()
        {
            if (filename.empty()) { std::cout << *this; std::cout.flush(); }
            else { std::ofstream file(filename); file << *this; }
            this->clear();
        }
    };

}
#endif // PRG_LOGGER
