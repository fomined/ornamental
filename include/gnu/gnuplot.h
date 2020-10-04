#ifndef PRG_GNUPLOT_H
#define PRG_GNUPLOT_H

#include <cstddef>
#include <string>
#include <vector>

#include "gnu/colors.h"
#include "math/curve.h"

namespace prg
{
    extern std::string TMP_DIR; // глобальная переменная, определяющая
    // директорию записи временных файлов

    /// пакет команд для выполнения gnuplot
    struct gnuplot_packet : public std::string
    {
        gnuplot_packet() = default;
        gnuplot_packet(const char *s) { this->append(s); }
        gnuplot_packet(const std::string &s) { this->append(s); }
        gnuplot_packet(const gnuplot_packet &) = default;
    
        bool empty() const { return std::string::empty(); }
    
        gnuplot_packet &operator<< (const char *s) { this->append(s); return *this; }
        gnuplot_packet &operator<< (const std::string &s) { this->append(s); return *this; }
        gnuplot_packet &operator<< (char s) { this->append(1, s); return *this; }
        gnuplot_packet &operator<< (const gnuplot_packet &p) { this->append(p); return *this; }
    
        gnuplot_packet operator()(const char *s) { gnuplot_packet g(*this); g.append(s); return g; }
        gnuplot_packet operator()(const std::string &s) { gnuplot_packet g(*this); g.append(s); return g; }
        gnuplot_packet operator()(char s) { gnuplot_packet g(*this); g.append(1, s); return g; }
    
        template <typename T> gnuplot_packet &operator<< (T s) { this->append(std::to_string(s)); return *this; }
        template <typename T> gnuplot_packet &operator()(T s) { this->append(std::to_string(s)); return *this; }
    };
    
    template <typename T>
    inline gnuplot_packet operator+(const gnuplot_packet &a, const T &b)
        { return gnuplot_packet(a) << b; }
    
    template <typename T>
    inline gnuplot_packet operator+(const T &a, const gnuplot_packet &b)
        { return gnuplot_packet(a) << b; }
    
    /// обертка интерпретатора gnuplot
    class gnuplot
    {
        FILE *pipe = nullptr;
    
    public:
    
        gnuplot() { pipe = popen ("gnuplot -p 2>/dev/null", "w"); } // избегаем предупреждений gnuplot
        ~gnuplot() { fflush(pipe); pclose(pipe); }
        gnuplot &operator<< (const gnuplot_packet &s) { fprintf(pipe, "%s\n", s.c_str()); return *this; }
    };
    
    /*!
     * \brief сбросить данные во временный файл для отрисовки их далее
     * \param i номер временого файла
     * \param p отрисовываемая кривая
     */
    template <Unsigned n>
    inline void gnuplot_prepare(Unsigned i, shared<curve<n>> p)
    {
        std::string TEMP_FILE = TMP_DIR + "/pattern.temp";
        FILE *temp = fopen((TEMP_FILE + std::to_string(i)).c_str(), "w");
        for (auto e : *p)
        {
            for (Unsigned k=0; k<n; k++) std::fprintf(temp, "%4.3lf ", e[k]);
            std::fprintf(temp, "\n");
        }
        fclose(temp);
    }

    /*!
     * \brief вставить данные в текущий график
     * \param pipe сценарий gnuplot
     * \param patterns данные по кривым
     * \param options опции отрисоки кривых
     */
    template <Unsigned n>
    inline void gnuplot_insert_plot(gnuplot &pipe,
        const std::vector<shared<curve<n>>> &patterns,
        const std::vector<gnuplot_packet> &options)
    {
        std::string TEMP_FILE = TMP_DIR + "/pattern.temp";

        gnuplot_packet gnu_plot("plot");
        Unsigned ns = patterns.size();
        for (Unsigned i=0; i<ns; i++)
        {
            if (i != 0) gnu_plot << ",";
            std::string tmpfile = std::string("'") + TEMP_FILE + std::to_string(i) + "' ";
            gnu_plot << tmpfile << options[i];
        }
        gnu_plot << ", 0 w l lc 'black' title '' "; // линия нуля вдоль x
        pipe << gnu_plot;
    }

    /*!
     * \brief стандартная рекламная тестовая диаграмма от gnuplot
     */
    inline void gnuplot_test(const std::string &output)
    {
        //---------------------------------------------------------------------------------
        //              рисуем диаграмму с помощью вызова gnuplot
        //---------------------------------------------------------------------------------
   
        gnuplot pipe; // инициализация gnuplot
        if (!output.empty())
        {
            pipe << std::string("set output '" + output + "'");
            pipe << "set term pdf enhanced";
        }
        pipe << "test";
    
        if (!output.empty())
        {
            printf("gnuplot_diagram '%s' is drawn\n", output.c_str());
            std::cout.flush();
        }
    }

} // namespace
#endif // PRG_GNUPLOT_H
