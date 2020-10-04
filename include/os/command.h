#ifndef OS_COMMAND_H
#define OS_COMMAND_H

/**
* \file
* \brief Команды, переадресующие operator() на функции выполнения
*/
#include <cstddef>
#include <stdlib.h>
#include <limits.h>
#include <cassert>
#include <string>
#include <map>
#include <queue>
#include <type_traits>
#include <fstream>
#include <iostream>
#include <ncurses.h> // нестандартное расширение для unix

#include "defs.h"
#include "container/string.h"
#include "os/os.h"
#include "math/maxmin.h"

//apropos -s 1 pdf
//dpkg --list | grep -i '[ (]pdf[ )].*viewer';
//apt install libncursesw5-dev

namespace prg
{

    enum { SUCCESS = 0, REINIT_NCURCES, REDRAW };

    /*!
     *\brief записать в неограниченный журнал сообщение (команду, ее результат)
     * \param extension расширение у названия файла
     * \note ".user" файл команд + ошибки выполнения (если есть, иначе пусто)
     * \note ".log" файл команд + результат выполнения
     */
    inline void filelog(const std::string &msg, const char *extension)
    {
        if (msg.empty()) return;

        extern std::string TMP_DIR;
        std::string filename = TMP_DIR + "/terminal" + extension;
        std::ofstream file(filename, std::ios_base::app);
        file << msg << "\n";
    }

    /*!
     *\brief записать в файл терминала сообщение (команду, ее результат)
     * \param extension расширение у названия файла
     * \note ".user" файл команд + ошибки выполнения (если есть, иначе пусто)
     * \note ".log" файл команд + результат выполнения
     */
    inline void terminal(const std::string &msg)
    {
        if (msg.empty()) return;

        extern std::string TMP_DIR;
        std::string filename = TMP_DIR + "/terminal.shot";
        std::ofstream file(filename);
        file << msg << "\n";
    }

    inline int cmd_msg(const std::string &msg)
    {
        terminal(msg); // вывод на экран
        filelog(msg, ".log"); // вывод в журнал
        return REDRAW;
    }

    inline int cmd_err(const std::string &msg)
    {
        cmd_msg(msg);
        return REDRAW;
    }

    inline int SYSTEM(const std::string &command)
    {
        if (command.empty()) return 0;
        filelog(command, ".log");

        extern std::string TMP_DIR;
        std::string txtfile = TMP_DIR + "/terminal.shot";
        std::string logfile = TMP_DIR + "/terminal.log";

        std::string cmd = command;
        cmd += " >" + txtfile + " 2>>" + logfile;
            // результат вызова системной команды пишется только на экран, но не в log файл
            // в log файл только ошибки

        auto r = system(&cmd[0]);
        if (r != EXIT_SUCCESS)
        {
            std::string msg = "Error : invalid command :\n";
            msg += command;
            terminal(msg);
            return REDRAW;
        }
        return REDRAW;
    }


    extern std::string TMP_DIR; // temporary files directory

    /// общий вид вызова функции, на которую настроен проект с пересылками команд
    typedef int (*PF)(const std::string &);

    /// список нужен для обеспечения безопасности (избежать случайного набора опасной команды)
    const std::set<std::string> system_commands
    {
        "pwd", "ls", "man",
        "apropos", "dpkg"
    };

    /*!
     * \brief команда выполнения
     */
    class Command
    {
        PF f_ = nullptr; // функция, которая вызывается командой (без функции - недействительна)
        std::string command = ""; // полная строка команды (с именем)

        /// связать с функций выполняющую команду из базы данных
        void assign(const std::map<std::string, PF> *addresses)
        {
            std::string name = this->name();
            f_ = addresses->count(name) ? addresses->at(name) : nullptr;
        }

        /// вызов системной команды (с необходимыми преобразованиями вывода под ncurses)
        int run_system(const std::string &command)
        {
            if (command.empty()) return 0;
            return SYSTEM(command);
        }

    public:

        Command(const std::string &s, const std::map<std::string, PF> *addresses) : command(trim(s))
            { assign(addresses); }

        /// вернуть имя команды для возможности связки с адресом функции по имени
        std::string name() const
        {
            if (command.empty()) return "";
            auto pos = command.find_first_of(' ');
            return std::string(command, 0, pos);
        }

        std::string param() const
        {
            assert( !command.empty() );
            auto pos = command.find_first_of(' ');

            if (pos == std::string::npos) return "";
            std::string p(command, pos + 1);

            return ltrim(p);
        }

        /// выполнение комманды
        int operator()() { return f_ ? f_(param()) : run_system(command); }

    };

    //-----------------------------------------------------------------------------------
    //                                  COMMAND FILE
    //-----------------------------------------------------------------------------------

    class cmd_Queue : public std::queue<Command>
    {
        const std::map<std::string, PF> *addresses; // таблица адресов функций
        WINDOW *mainwin = nullptr; // main window of command terminal
        WINDOW *cmdwin = nullptr; // part of main window (commands)
        WINDOW *outwin = nullptr; // part of main window (results)

        void open_command_terminal();
        void close_command_terminal();

        std::string get_command(const accumulator<std::string> &history, WINDOW *win);
        std::string get_command_name(const std::string &command);

        int create_tmp_directory();
        int delete_tmp_directory();

        // переинсталлируем ncurses с новой локалью
        void reinit_ncurces()
        {
            close_command_terminal();
            open_command_terminal();
        }

        // вывод данных, полученных после выполнения команды, на терминал
        void output();


    public:

        cmd_Queue(const std::map<std::string, PF> *addresses, const std::string &filename="")
            : addresses(addresses)
        {
            create_tmp_directory(); // создать директорию для временных файлов пользователя
            load(filename); // load init file
        }

        ~cmd_Queue()
        {
            close_command_terminal();
            // delete_tmp_directory();
        }

        /// load init commands from default directory [~/.ornamental]
        void load(const std::string &filename);

        /// получать команды из консоли и выполнять их пока не появится exit()
        void cycle();
    };

    inline std::string cmd_Queue::get_command(const accumulator<std::string> &history, WINDOW *win)
    {
        Signed n = history.size();

        std::string new_command;

        Unsigned pos = new_command.size();
        Signed row = -1; // позиция новой команды

        wclear(win);
        wrefresh(win);

        int ch = 0;
        while (ch != '\n')
        {
            ch = wgetch(win);
            switch (ch)
            {
            case KEY_UP     :
                if (row < n - 1) ++row;
                break;

            case KEY_DOWN   :
                if (row > -1) --row; // возможность npos
                break;

            case KEY_LEFT   :
            case KEY_RIGHT  :
                break;

            case KEY_BACKSPACE  :
                if (row >= 0) { new_command = history[row]; row = -1; pos = new_command.size(); }
                if (pos != 0) new_command.erase(--pos, 1);
                break;

            default:
                if (isascii(ch) && isprint(ch))
                {
                    if (row >= 0) { new_command = history[row]; row = -1; pos = new_command.size(); }
                    new_command.insert(pos++, 1, (char)ch);
                }
            }

            wclear(win);
            wmove(win, 0, 0);
            if (row >= 0) waddstr(win, history[row].c_str());
            else waddstr(win, new_command.c_str());

            wrefresh(win);
        }

        if (row >= 0) return history[row];
        return new_command;

    }

    /// извлечь имя команды
    inline std::string cmd_Queue::get_command_name(const std::string &command)
    {
        if (command.empty()) return "";
        auto pos = command.find_first_of(' ');
        return std::string(command, 0, pos);
    }

    inline int cmd_Queue::create_tmp_directory()
    {
        extern std::string HOME_DIR;

        HOME_DIR = std::string(getenv("HOME"));
        TMP_DIR = HOME_DIR + "/.ornamental";

        if (1 != dir_exists(TMP_DIR.c_str()))
        {
            std::string cmd = "mkdir " + TMP_DIR;
            system(cmd.c_str()); // создать директорию временных файлов пользователя
        }

        std::string filename = TMP_DIR + "/terminal.log";
        std::string command("printf '' ");
        command += " >" + filename;
        system(&command[0]); // создать пустой log файл

        filename = TMP_DIR + "/terminal.user";
        command = "printf '' ";
        command += " >" + filename;
        system(&command[0]); // создать пустой cmd файл

        command = "rm " + TMP_DIR + "/terminal.svg " + TMP_DIR + "/terminal.pdf 2>/dev/null";
        system(&command[0]); // удалить данные с прежнего сеанса

        return 0;
    }

    inline int cmd_Queue::delete_tmp_directory()
    {
        std::string cmd = "rm -r " + TMP_DIR;
        return system(cmd.c_str());
    }

    inline void cmd_Queue::output()
    {
        int ymax = getmaxy(outwin);
        wclear(outwin);

        // вывод данных, полученных после выполнения команды, на терминал
        char buffer[1024];
        std::ifstream file(TMP_DIR + "/terminal.shot");

        int y=0, x=0;
        file.getline(buffer, 1024);
        while (!file.eof())
        {
            std::string s(buffer);

            wmove(outwin, y++, x);
            waddstr(outwin, s.c_str());

            if (y == ymax)
            {
                wrefresh(outwin);
                wgetch(outwin);
                wclear(outwin);
                wrefresh(outwin);
                y = 0;
            }

            file.getline(buffer, 1024);
        }
        wrefresh(outwin);
    }

    inline void cmd_Queue::load(const std::string &filename)
    {
        const Unsigned LEN = 1024;
        char buffer[LEN];

        std::ifstream file(filename);
        if (!file)
        {
            std::string msg = "cmd_Queue::load() : \'" + filename;
            msg += "\' file opening failed\n";
            filelog("cmd_Queue::load()", ".log");
            filelog(msg, ".log");
        }
        else
        {
            file.getline(buffer, LEN);
            while (!file.eof())
            {
                std::string command(buffer);
                auto pos = command.find_first_of('#');
                if (pos != std::string::npos) command.erase(pos);

                filelog(command, ".user");
                filelog(command, ".log");
                std::string name = get_command_name(command);
                if (name == "exit" || name == "quit") exit(0);

                int r = Command(command, addresses)(); // создание и выполнение команды
                if (r == REINIT_NCURCES) reinit_ncurces();

                file.getline(buffer, 1024);
            }
        }
        open_command_terminal();
    }

    inline void cmd_Queue::open_command_terminal()
    {
        close_command_terminal();

        std::string command = "";

        extern std::string LOCALE;
        std::string locale = LOCALE; // to view
        std::setlocale(LC_ALL, LOCALE.c_str());

        mainwin = initscr(); // determines the terminal type, including its size and features,
            // and sets up the curses environment

        if (mainwin == nullptr)
        {
            filelog("cmd_Queue::open_command_terminal()", ".log");
            filelog("Error initialising ncurses.", ".log");
            exit(EXIT_FAILURE);
        }

        cbreak(); // disables line buffering
        // noecho(); // not to echo the input back to the screen
        clear(); // clears the screen
        keypad(stdscr, TRUE); // to enable support for function-key mapping

        std::string prompt = "Enter command: ";
        int startx = prompt.length();

        int height, width;
        getmaxyx(stdscr, height, width); // get the size of the screen

        cmdwin = newwin(1, width-startx, 1, startx);
        outwin = newwin(height - 4, width, 2, 0);
//        infowin = newwin(1, width, height-2, 0);

        keypad(cmdwin, TRUE);
        keypad(outwin, TRUE);
//        keypad(infowin, TRUE);

        mvaddstr(0, 0, "if the terminal does not respond press \'Ctrl+C\'");
        mvaddstr(1, 0, prompt.c_str()); // prompt is placed out of all windows
        refresh();
    }

    inline void cmd_Queue::close_command_terminal()
    {
        if (cmdwin) delwin(cmdwin);
        cmdwin = nullptr;

        if (outwin) delwin(outwin);
        outwin = nullptr;

//        if (infowin) delwin(infowin);
//        infowin = nullptr;

        // exit the curses environment and return the terminal screen to normal control
        if (mainwin) endwin();
        mainwin = nullptr;
    }

    inline void cmd_Queue::cycle()
    {
        std::string command = "";

        History<std::string> history(15);
        while (true)
        {
            std::string command = get_command(history, cmdwin);
            filelog(command, ".user");
            filelog(command, ".log");

            std::string name = get_command_name(command);
            if (name == "exit" || name == "quit")
            {
                wclear(outwin);
                wrefresh(outwin);
                break;
            }
            if (!command.empty()) history.push(command);

//            wmove(infowin, 0, 0); // перейдем в окно ошибок для возможной их печати
//            wrefresh(infowin);

            int r = Command(command, addresses)(); // создание и выполнение команды
            if (r == REINIT_NCURCES) reinit_ncurces();
            if (r == REDRAW) output();

//            wrefresh(infowin);
        }
    }

}
#endif // OS_COMMAND_H
