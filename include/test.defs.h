#ifndef TEST_DEFS_H
#define TEST_DEFS_H

#include <cstddef>
#include <numeric>
#include <array>
#include <vector>
#include <iostream>

/*!
 * \file файл стандартных определений типов и функций, используемых для отладки и тестирования
 */

#define TEST_EQ(a, b) \
    { \
        bool expr = fabs(a - b) < 0.001 * std::min(fabs(a), fabs(b)); \
        if (!expr) std::cout << FILE_LINE << " failed ... \n"; \
    }

#define TEST_LE(a, b) \
    { \
        if (a > b) std::cout << FILE_LINE << " failed ... \n"; \
    }

template <typename T>
inline void print_compare(const char *m, T t1, T t2)
{ std::cout << m << " : " << t1 << " vs " << t2 << "\n"; std::cout.flush(); }


/// печатать разделительную линию с заданным символом в output
inline void print_separator(size_t n=80, char ch='-')
{ printf("%s\n", std::string(n, ch).c_str()); std::cout.flush(); }


#endif // TEST_DEFS_H
