#ifndef OS_H
#define OS_H

#include <fstream>
#include <string>
#include <sys/stat.h>

#ifdef BOOST_USED
    #include <boost/filesystem/path.hpp>
    #include <boost/filesystem/convenience.hpp>
#endif

namespace prg
{

    constexpr size_t _1K   = 1024;
    constexpr size_t _1M   = 1024 * _1K;
    constexpr size_t _1G   = 1024 * _1M;

    /**
     * @brief превратить имя в имя директории, то есть добавить "/" в конце при необходимости
     * @note гарантирует, что символ "/" будет в конце и будет единственным
     * @param name строка
     * @return строка с "/" в конце
     */
    inline std::string make_dirname(const std::string &name)
    {
        if (name.empty()) return std::string("./");
        if (name[name.size() - 1] == '/') return name;
        return name + std::string("/");
    }

    /*!
     * \brief checks to see if a directory exists.
     * \note This method only checks the existence of the full path AND if path leaf is a dir.
     * \return  1 if dir exists AND is a dir,
     *     0 if dir does not exist OR exists but not a dir,
     *    -1 if an error occurred (errno is also set)
     */
   inline int dir_exists(const char* path)
   {
       struct stat info;
       int result = stat(path, &info);
       if (result != 0)
       {
           if (errno == ENOENT)  { return 0; } // not exist
           if (errno == ENOTDIR) { return 0; } // not a dir
           return -1; // error
       }
       return (info.st_mode & S_IFDIR) ? 1 : 0;
   }

#ifdef BOOST_USED

   inline std::string extension(const std::string &s)
   { return boost::filesystem::extension(boost::filesystem::path(s)); }

    /// получить имя файла без расширения
    inline std::string basename(const std::string &s)
    { return boost::filesystem::basename(boost::filesystem::path(s)); }

    /// вычленить директорию из имени файла
    inline std::string parent_path(const std::string &s)
    { return boost::filesystem::path(s).parent_path().native(); }

    /// вычленить имя файла из пути
    inline std::string filename(const std::string &s)
    { return boost::filesystem::path(s).filename().native(); }

#else

    /// get extension of file (.*) only for unix
    inline std::string extension(const std::string &s)
    {
        auto npos = std::string::npos;

        auto pos2 = s.find_last_of('.');
        if (pos2 == npos) return "";

        auto pos1 = s.find_last_of('/');
        if (pos1 == npos || pos1 < pos2) return std::string(s, pos2);

        return "";
    }

    /// получить имя файла без расширения
    inline std::string basename(const std::string &s)
    {
        auto npos = std::string::npos;

        auto pos1 = s.find_last_of('/');
        if (pos1 == npos) pos1 = 0;

        auto pos2 = s.find_last_of('.');
        if (pos2 == npos) pos2 = s.length();

        return pos1 < pos2 ? std::string(s, pos1, pos2 - pos1) : std::string(s, pos1);
    }

#endif

}
#endif // OS_H
