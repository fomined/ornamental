#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>

#include "prg.os.h"
#include "prg.string.h"
#include "flower.structure.h"

bool Flowerstruct::load_file(const std::string &filename, bool)
{
    if (filename.empty()) return false;
    species_name = prg::basename(filename); // сохраним имя вида

    char buffer[BUF_LEN], buf[BUF_LEN];

    std::ifstream file(filename + "-axe.csv");
    if (!file) no_file_exit(filename);

    file.getline(buffer, BUF_LEN); // пропуск строки с именами колонок

    unsigned prev_year = 0;
    file.getline(buffer, BUF_LEN);
    while (!file.eof())
    {
        std::stringstream ss(buffer);

        ss.getline(buf, BUF_LEN, ',');
        char *p = prg::trim(buf, remove_symbols);
        unsigned year = (strlen(p) == 0) ? prev_year : unsigned(atoi(p));
        prev_year = year;

        ss.getline(buf, BUF_LEN, ',');
        p = prg::trim(buf, remove_symbols);
        unsigned sprout = unsigned(atoi(p));

        ss.getline(buf, BUF_LEN, ',');
        p = prg::trim(buf, remove_symbols);
        unsigned n = unsigned(atoi(p));

        ss.getline(buf, BUF_LEN, ',');
        p = prg::trim(buf, remove_symbols);
        unsigned n1 = unsigned(atoi(p));

        push_back({year, sprout, n1, n - n1});

        file.getline(buffer, BUF_LEN);
    }
    return true;
}

curve_t Flowerstruct::count_select(const std::vector<unsigned> &years, unsigned axis) const
{
    std::map<unsigned, unsigned> num;

    for (auto e : years)
    {
        for (auto i : *this)
        {
            if (i.year != e) continue;
            auto n = (axis == 1) ? i.n1 : i.n2;
            if (num.count(n) == 0) num[n] = 0;
            num[n]++;
        }
    }

    curve_t data;
    for (auto e : num) data.push_back({double(e.first), double(e.second)});
    return data;
}

curve_t Flowerstruct::count_select(const std::vector<unsigned> &years) const
{
    std::map<unsigned, unsigned> num;

    for (auto e : years)
    {
        for (auto i : *this)
        {
            if (i.year != e) continue;
            auto n = i.n1 + i.n2;
            if (num.count(n) == 0) num[n] = 0;
            num[n]++;
        }
    }

    curve_t data;
    for (auto e : num) data.push_back({double(e.first), double(e.second)});
    return data;
}

curve_t Flowerstruct::rate_select(const std::vector<unsigned> &years) const
{
    curve_t data;

    for (auto e : years)
    {
        for (auto i : *this)
        {
            if (i.year != e) continue;
            data.push_back({double(i.n1), double(i.n2)});
        }
    }
    std::sort(data.begin(), data.end(),
        [](const point_t &a, const point_t &b)
        { return a[0] < b[0]; }
    );

    return data;
}


/////////////////////////////////////////////////////////////////////////////////////////
//curve_t Flowerstruct::count_select(const std::vector<unsigned> &years) const
//{
//    std::map<double, unsigned> num;

//    for (auto e : years)
//    {
//        for (auto i : *this)
//        {
//            if (i.year != e) continue;
//            auto n = double(i.n2) / i.n1;
//            if (num.count(n) == 0) num[n] = 0;
//            num[n]++;
//        }
//    }

//    curve_t data;
//    for (auto e : num) data.push_back({double(e.first), double(e.second)});
//    return data;
//}

