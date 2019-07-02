#ifndef OS__H
#define OS__H

#include <fstream>
#include <string>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>

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

    /// получить расширение файла (.*)
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

    /// итератор по директории
    class directory_iterator : public boost::filesystem::directory_iterator
    {
        typedef boost::filesystem::directory_iterator base;
    public:

        /// iterator begin()
        directory_iterator(const std::string &dirname) : base(boost::filesystem::path(dirname)) {}

        /// iterator end()
        directory_iterator() : base() {}

        /// iterator of extraction of file name
        std::string operator*() const { return ((*this)->path()).filename().generic_string(); }

        /// iterator of moving to next file name (directories will bw skipped)
        directory_iterator &operator++()
        {
            do { base::operator++(); }
            while ( (*this) != base() && boost::filesystem::is_directory(**this) );
            return *this;
        }

        /// оператор сравнения
        bool operator==(const directory_iterator &it) const
        {
            return (const boost::filesystem::directory_iterator&)(*this)
                == (const boost::filesystem::directory_iterator&)it;
        }

        /// оператор сравнения
        bool operator!=(const directory_iterator &it) const { return !(*this == it); }
    };

    /// пара итераторов (range) по директории
    class directory_range : std::pair<directory_iterator, directory_iterator>
    {
        typedef std::pair<directory_iterator, directory_iterator> base;

    public:

        directory_range(const std::string &dirname)
        : base(directory_iterator(dirname), directory_iterator()) {}

        /// iterator of extraction of file name
        std::string front() const { return *base::first; }

        /// iterator of moving to next file name (directories will be skipped)
        directory_range &pop_front() { ++base::first; return *this; }

        /// проверка завершения области итерирования
        bool empty() const { return base::first == base::second; }
    };


}
#endif // OS__H
