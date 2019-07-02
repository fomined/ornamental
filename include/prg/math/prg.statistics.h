#ifndef PRG_STATISTICS_H
#define PRG_STATISTICS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

#include "prg.defs.h"

namespace  prg
{
    /// гистограмма
    class Histogram : protected std::vector<Unsigned>
    {
        typedef std::vector<Unsigned> base_t;
        double a, b, dx; // границы и ширина корзины
    
        flagged_value<Unsigned> index(double v) const
        { 
            if (a <= v && v < b)
            {
                Unsigned k = Unsigned((v - a) / dx); 
                assert(k < size()); 
                return {k, true}; 
            }
            return {0, false};
        }
    
        double value(Unsigned k) const { return k * dx; }
    
    public :
    
        Histogram(double a, double b, double dx) 
        : base_t(Unsigned((b - a) / dx) + 1, 0), a(a), b(b), dx(dx) {}

        Histogram(double a, double b, Unsigned n) 
        : base_t(n, 0), a(a), b(b), dx((b - a)/n) {}
    
        /// вставка только тех значений, которые находятся в диапазоне [a, b)
        bool insert(double v) { auto f = index(v); if (f.second) at(f.first)++; return f.second; }
    
        void print()
        {
            std::printf("--------- histogram start ---------\n");
            for (Unsigned i=0; i<size(); i++)
            {
                std::printf("%lf %ld\n", value(i), at(i));
            }
            std::printf("--------- histogram end ---------\n");
            FLUSH;
        }
    };

    /// получение среднего для массива
    inline double average(const std::vector<double> &v)
    {
        size_t n = v.size();
        double sum = 0; for (size_t i=0; i<n; i++) sum += v[i];
        return sum / n;
    }

    /// получение дисперсии для массива
    /// !делать по умалчиванию mean=0 нельзя, так как провоцирует ошибки забывания
    inline double variance(const std::vector<double> &v, double mean)
    {
        size_t n = v.size();
        double sum = 0; for (size_t i=0; i<n; i++) sum += sqr(v[i] - mean);
        return sum / n;
    }

} // namespace prg
#endif // PRG_STATISTICS_H
