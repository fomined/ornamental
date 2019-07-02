#include <string>
#include <vector>
#include <set>
#include <array>
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <cctype>
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>

#include "prg.random.h"
#include "prg.string.h"
#include "prg.os.h"
#include "defs.h"
#include "flower.dome.h"

void test_file(const std::string &filename)
{
    Flowerdome dome;
    dome.load_file(filename, {}, true);
}

void test_directory(const std::string &dirname)
{
    using namespace boost::filesystem;
    for (auto &entry : boost::make_iterator_range(directory_iterator(dirname), {}))
    {
        std::string filename = entry.path().filename().native();
        std::string ext = prg::extension(filename);
        if (ext != ".csv") continue;

        printf("\n################ testing of file \'%s\' ################\n", filename.c_str());
        test_file(dirname + filename);
    }
}
