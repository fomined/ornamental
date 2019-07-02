#ifndef RECORD_TIME_H
#define RECORD_TIME_H

#include <string>
#include <iostream>
#include <cctype>

#include "prg.string.h"
#include "defs.h"

/// частичный контроль того, является ли поле представлением даты
inline bool is_date(const std::string &s)
{
    std::string::size_type pos = s.find("-");
    if (pos == std::string::npos || pos == s.length()-1) return false;

    for (size_t i=0; i<pos; i++)
        if (!std::isdigit(static_cast<unsigned char>(s[i]))) return false;

    if (strncmp(&s[pos+1], "Jun", 3) == 0) return true;
    if (strncmp(&s[pos+1], "Jul", 3) == 0) return true;
    if (strncmp(&s[pos+1], "Aug", 3) == 0) return true;
    if (strncmp(&s[pos+1], "Sep", 3) == 0) return true;

    return false;
}

/// облегченная структура хранения времени (день лета, т.е. отчет с 1 июня)
class Timestamp
{
    int tm;

public:

    /// получить время из года и строки формата "1-Jun"
    Timestamp(const char *s)
    {
        char buffer[128]; std::strcpy(buffer, s);
        const char *p = prg::trim(buffer, remove_symbols);
        tm = atoi(p); // день месяца

        const char *q = strpbrk(p, "-"); q++;
        if (strncmp(q, "Jun", 3) != 0) tm += 30; else return;
        if (strncmp(q, "Jul", 3) != 0) tm += 31; else return;
        if (strncmp(q, "Aug", 3) != 0) tm += 31; else return;
        if (strncmp(q, "Sep", 3) != 0) EXIT;
    }

    /// устанавливается объект с известным днем лета (либо 0)
    Timestamp() : tm(0) {}
    Timestamp(int s) : tm(s) {}
    Timestamp(unsigned s) : tm(int(s)) {}

    /// устанавливается объект с известным днем лета (либо 0)
    Timestamp(double s) : tm(int(round(s))) {}

    /// оператор сдвига на следующую дату
    Timestamp &operator++() { tm += 1; return *this; }

    /// конвертация для сравнения по дням
    int to_int() const { return tm; }
    unsigned to_unsigned() const { return unsigned(tm); }
    double to_double() const { return double(tm); }

    /// действительная дата
    bool valid() const { return tm > 0; }

    /// конвертация в строку "1-Jul"
    std::string to_string() const
    {
        int t = tm;
        if (t <= 0) return std::string("invalid");
        if (t <= 30) return std::to_string(t) + std::string("-Jun");
        t -= 30; if (t <= 31) return std::to_string(t) + std::string("-Jul");
        t -= 31; if (t <= 31) return std::to_string(t) + std::string("-Aug");
        t -= 31; if (t <= 30) return std::to_string(t) + std::string("-Sep");
        t -= 30; if (t <= 31) return std::to_string(t) + std::string("-Oct");
        return std::string("invalid");
    }

    Timestamp &operator+=(int n) { tm += n; return *this; }
    Timestamp &operator-=(int n) { tm -= n; return *this; }
};

/// вывод на печать метки времени
inline std::ostream &operator <<(std::ostream & os, const Timestamp &s)
{ return os << s.to_string(); }

inline Timestamp operator+(const Timestamp &t1, int t2) { return t1.to_int() + t2; }
inline Timestamp operator-(const Timestamp &t1, int t2) { return t1.to_int() - t2; }

/// протяженность временного интервала
inline int operator-(const Timestamp &t1, const Timestamp &t2)
{ return t1.to_int() - t2.to_int(); }

inline bool operator==(const Timestamp &t1, const Timestamp &t2)
{ return t1.to_int() == t2.to_int(); }

inline bool operator!=(const Timestamp &t1, const Timestamp &t2)
{ return t1.to_int() != t2.to_int(); }

inline bool operator<(const Timestamp &t1, const Timestamp &t2)
{ return t1.to_int() < t2.to_int(); }

inline bool operator<=(const Timestamp &t1, const Timestamp &t2)
{ return t1.to_int() <= t2.to_int(); }

inline bool operator>(const Timestamp &t1, const Timestamp &t2)
{ return t1.to_int() > t2.to_int(); }

inline bool operator>=(const Timestamp &t1, const Timestamp &t2)
{ return t1.to_int() >= t2.to_int(); }

#endif // RECORD_TIME_H
