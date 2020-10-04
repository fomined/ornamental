#ifndef PRG_STRING
#define PRG_STRING

/**
* @file Отображение памяти в файл и чтение памяти из файла.
*/

#include <cstdio>
#include <cstring>
#include <sstream>
#include <string>
#include <memory>
#include <cctype>
#include "defs.h"

namespace prg
{

    /// литеральная строка со всеми пустыми символами
    const char * const EMPTY_SYMBOLS = " \t\n\r\f\v";

    /// является ли заданный символ одним из перечисленных
    inline bool is_oneof(char s, const char* t)
    {
        while (*t && (s != *t)) ++t;
        return *t != 0;
    }

    /// удалить ведущие пробелы в строке
    inline char *ltrim(char *s, const char* t=EMPTY_SYMBOLS)
    {
        while (*s && is_oneof(*s, t) ) ++s;
        return s;
    }

    /// удалить заключительные пробелы в строке
    inline char *rtrim(char *s, const char* t=EMPTY_SYMBOLS)
    {
        char *q = s + strlen(s);
        while (q != s && is_oneof(*(q-1), t)) --q;
        *q = 0; // обрезаем строку
        return s;
    }

    /// удалить ведущие и заключительные пробелы в строке
    inline char *trim(char *s, const char* t=EMPTY_SYMBOLS) { return ltrim(rtrim(s, t), t); }

    /// определение того, что строка не содержит значимых символов
    inline bool empty(char *s, const char* t=EMPTY_SYMBOLS)
    {
        char *p = ltrim(s, t);
        return unsigned(p - s) == strlen(s);
    }

    /// удалить ведущие пробелы в строке
    inline std::string ltrim(const std::string &s, const char* t=EMPTY_SYMBOLS)
    {
        std::string ss(s);
        ss.erase(0, ss.find_first_not_of(t));
        return ss;
    }

    /// удалить заключительные пробелы в строке
    inline std::string rtrim(const std::string &s, const char* t=EMPTY_SYMBOLS)
    {
        std::string ss(s);
        ss.erase(ss.find_last_not_of(t) + 1);
        return ss;
    }

    /// удалить ведущие и заключительные пробелы в строке
    inline std::string trim(const std::string &s, const char* t=EMPTY_SYMBOLS)
    {
        std::string ss(s);
        ss = ltrim(rtrim(ss, t), t);
        return ss;
    }

    /// форматный вывод в строку, аналог sprintf(char *, const char *, ..)
    inline std::string to_string(const char *format)
    {
        std::string msg;
        while (*format)
        {
            if (*format == '%' && *(format+1) != '%')
                throw std::runtime_error("invalid format string: missing arguments");
            msg += *format++;
        }
        return msg;
    }

    /// форматный вывод в строку, аналог sprintf(char *, const char *, ..)
    template<typename T, typename... Args>
    inline std::string to_string(const char *format, T value, Args... args)
    {
        std::string msg;
        while (*format)
        {
            if (*format == '%' && *(format+1) != '%')
            {
                char argument_format[64];

                const char *p = format; int i = 0;
                argument_format[i++] = *p; p++;
                while (!isspace(*p) && *p!=0 && *p!='%')
                {
                    argument_format[i++] = *p;
                    ++p;
                }
                argument_format[i] = '\0';

                char argument_value[1024];
                std::sprintf(argument_value, argument_format, value);
                msg.append(argument_value);

                format = p;
                msg += to_string(format, args...);
                return msg;
            }
            msg += *format++;
        }
        throw std::logic_error("extra arguments provided to printf");
    }

    inline std::vector<std::string> split(const std::string &line, char delim)
    {
        std::vector<std::string> out;
        size_t n = line.length();

        std::istringstream iss(line); char buf[n + 1];
        for (; iss.getline(buf, static_cast<std::streamsize>(n+1), delim); )
            // n + 1, для считывания конца строки (иначе теряется первый символ)
        {
            std::string v = (buf[0] != '\0' && buf[0] != delim) ? buf : "";
            if (!v.empty()) out.push_back(std::string(buf));
        }
        return out;
    }

    inline std::vector<std::string> split(const std::string &s, const std::string &delims=std::string(EMPTY_SYMBOLS))
    {
        std::vector<std::string> out;

        Unsigned p1 = 0, p2 = 0;
        while (true)
        {
            p1 = s.find_first_not_of(delims, p2);
            if (p1 == std::string::npos) return out;

            p2 = s.find_first_of(delims, p1);
            out.push_back(std::string(s, p1, p2 == std::string::npos ? p2 : p2 - p1));
        }
        return out;
    }

//    template <typename T>
//    inline std::vector<T> split(const std::string &line, T empty, char delim=',')
//    {
//        std::vector<T> r;
//        size_t n = line.length();

//        std::istringstream iss(line); std::string buf(n + 1, ' ');
//        for (; iss.getline(&buf.at(0), static_cast<std::streamsize>(n+1), delim); )
//            // n + 1, для считывания конца строки (иначе теряется первый символ)
//        {
//            T v = (buf.at(0) != '\0') ? boost::lexical_cast<T>(&buf.at(0)) : empty;
//            r.push_back(v);
//        }
//        if (line.back() == delim) r.push_back(empty);

//        return r;
//    }

    /*!
     * \brief конвертировать число в строку с заданным числом десятичных разрядов
     * \note функция работает для любого типа (кроме строки с точкой внутри)
     * \param W ширина поля вывода
     * \param w число разрядов после запятой
     */
    template <typename T>
    inline std::string to_string(T a, unsigned W=6, unsigned w=3) 
    {
        auto s = std::to_string(a);
        
        // отрезаем конец строки после нужного числа разрядов 
        auto pos = s.find_first_of('.');
        if (pos != std::string::npos) { pos += w; s.resize(pos + 1, '0'); }
        
        // вставим незначащие пробелы впереди
        if (s.size() < W) s.insert(s.begin(), W - s.size(), ' ');
        
        return s;
    }

}
#endif // PRG_STRING
