#ifndef ORN_TIME_H
#define ORN_TIME_H

#include <string>
#include <iostream>
#include <cctype>

#include "container/vector.h"
#include "orn/orn.defs.h"

namespace prg
{
    const double month_days[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    const char * const month_names[12] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};

    /// частичный контроль того, является ли поле представлением даты
    inline bool is_date(const std::string &p)
    {
        std::string s = trim(p, remove_symbols);
        std::string::size_type pos = s.find("-");
        if (pos == std::string::npos) return false;

        for (Unsigned i=0; i<pos; i++)
            if (!std::isdigit(static_cast<unsigned char>(s[i]))) return false;
            // обязятельно преобразование в unsigned char (см. https://en.cppreference.com/)

        for (Unsigned i=0; i<12; i++)
            if (strncmp(&s[pos+1], month_names[i], 3) == 0) return true;

        return false;
    }

    inline Unsigned year_day(const std::string &p)
    {
        if (!is_date(p)) return npos;

        std::string s = trim(p, remove_symbols);
        std::string::size_type pos = s.find("-");

        double month_day = std::atof(&s[0]); // день месяца

        Unsigned month = npos;
        for (Unsigned i=0; i<12; i++)
            if (strncmp(&s[pos+1], month_names[i], 3) == 0) { month = i; break; }

        Unsigned sum = 0;
        for (Unsigned m=0; m<month; m++) sum += month_days[m];

        return sum + month_day;
    }

    const double invalid_time = -366.;

    /// облегченная структура хранения времени (день лета, т.е. отчет с 1 июня)
    class Time
    {
        double tm; // чтобы можно было оптимизировать

    public:
        typedef double difference_type;

        /// устанавливается объект с известным днем лета (либо 0)
        constexpr Time(double s=invalid_time) : tm(s) {}

        /// получить время из года и строки формата "1-Jun"
        explicit Time(const char *s) : tm(year_day(s)) {}

        /// действительная дата
        bool valid() const { return tm >= 0; }

        /// оператор сдвига на следующую дату
        Time &operator++() { tm += 1; return *this; }

        /// конвертация для сравнения по дням
        double to_double() const { return tm; }
        Signed to_signed() const { return (Signed)std::floor(tm); }

        Unsigned to_month_day() const
        {
            double t = tm;
            Unsigned i=0; while (t > month_days[i] && i<12) { t -= month_days[i]; i++; }
            if (i == 12) return 0;
            return t;
        }
        Unsigned to_month() const
        {
            double t = tm;
            Unsigned i=0; while (t > month_days[i] && i<12) { t -= month_days[i]; i++; }
            if (i == 12) return 0;
            return i+1;
        }

        /// конвертация в форматированную строку "1-Jul"
        std::string to_short_string() const
        {
            double t = tm;
            Unsigned i=0; while (t > month_days[i] && i<12) { t -= month_days[i]; i++; }
            if (i == 12) return "  ?  ";

            char buffer[128];
            std::sprintf(buffer, "%2ld.%ld", Unsigned(t), i+1);
            return buffer;
        }

        /// конвертация в форматированную строку "1-Jul"
        std::string to_string() const
        {
            double t = tm;
            Unsigned i=0; while (t > month_days[i] && i<12) { t -= month_days[i]; i++; }
            if (i == 12) return "   ?  ";

            char buffer[128];
            std::sprintf(buffer, "%2ld-%3s", Unsigned(t), month_names[i]);
            return buffer;
        }

        Time &operator+=(difference_type n) { tm += n; return *this; }
        Time &operator-=(difference_type n) { tm -= n; return *this; }
    };

    /// вывод в строку метки времени
    inline std::string to_string(const Time &s) { return s.to_string(); }

    /// вывод на печать метки времени
    inline std::ostream &operator <<(std::ostream & os, const Time &s)
    { return os << s.to_string(); }

    inline Time operator+(const Time &t1, Time::difference_type t2) { Time t(t1); return t += t2; }
    inline Time operator-(const Time &t1, Time::difference_type t2) { Time t(t1); return t -= t2; }

    /// протяженность временного интервала
    inline Time::difference_type operator-(const Time &t1, const Time &t2)
    { return t1.to_double() - t2.to_double(); }

    #define OPER(op) inline bool operator op(const Time &t1, const Time &t2) \
        { return t1.to_signed() op t2.to_signed(); }
        OPER(==)
        OPER(!=)
        OPER(<)
        OPER(<=)
        OPER(>)
        OPER(>=)
    #undef OPER

    inline constexpr Time::difference_type minimal_difference(Time) { return 1; }

    /// для поиска максимальных/минимальных дат
    template <> inline Time inf<Time>() { return Time(366.); }

}
#endif // RECORD_TIME_H
