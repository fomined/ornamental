#ifndef ORN_COMMAND_H
#define ORN_COMMAND_H

#include <cstddef>
#include <numeric>
#include <array>
#include <type_traits>
#include <sys/stat.h>

#include "defs.h"
#include "os/os.h"
#include "orn.defs.h"
#include "os/command.h"
#include "orn/flowerdome.h"
#include "orn/sproutmodel.h"

namespace prg
{
    // !!! команды не могут различаться по регистру!

    inline int cmd_help(const std::string &)
    {
        std::string msg;
        msg += "\n";
        msg += "help|h                                      help\n";
        msg += "environment|env [<var> <value>]             list of environment variables\n";
        msg += "window|w                                    open window to view the plots\n";
        msg += "------------------------------------------------------------------------------------\n";
        msg += "files|f                                     view status of the data files\n";
        msg += "files|f <fileID>                            view data file info\n";
        msg += "files|f <fileID> R|W|RW                     set status of the data file\n";
        msg += "fload|f< <filename> [<fileID>]              load sprout data file\n";
        msg += "fsave|f> <fileID> <filename>                save the file <ID>\n";
        msg += "fadd|f+ <fileID>                            create new data file\n";
        msg += "fdel|f- <fileID>                            delete the file <ID>\n";
        msg += "festimate|f^ MLM|random <fileID>            estimate of missing data\n";
        msg += "------------------------------------------------------------------------------------\n";
        msg += "models|m                                    view status of all the models\n";
        msg += "models|m <modelID>                          set the active model\n";
        msg += "mload|m< <filename> [<modelID>]             [reserved]\n";
        msg += "msave|m> <modelID> <filename>               [reserved]\n";
        msg += "madd|m+ <fileID:sproutID> [<modelID>]       build a model of the sprout\n";
        msg += "mdel|m- <modelID>                           delete the model modelID\n";
        msg += "------------------------------------------------------------------------------------\n";
        msg += "plots|p                                     show IDs of the all plots\n";
        msg += "plot|p <plotID>                             show the plots\n";
        msg += "padd|p+ <plotID> [width value]              add the plot with smoooth width\n";
        msg += "pdel|p-                                     delete the current plot\n";
        msg += "pdel|p- <plotID>                            delete the plot <ID>\n";
        msg += "psave|p> <filename>                         save the current plot (as pdf)\n";
        msg += "------------------------------------------------------------------------------------\n";
        msg += "scheme <modelID> [<date|day=0>]             show the scheme of inflorescence\n";
        msg += "scheme <fileID>:<sproutID> [<date|day=0>]   show the scheme of inflorescence\n";
        msg += "video <modelID> <date1:date2|day1:day2>     video of development of inflorescence\n";
        msg += "test|t <plotID> [nsamples]                  statistics and normality test for <plotID>\n";
        msg += "generate|g <modelID> <N>[<n1>:<n2>]         generate N model sprouts\n";
        msg += "\n";
        msg += "quit|exit                                   exit the program\n";

        return cmd_msg(msg);
    }

    inline bool test_format(const std::string &s, Unsigned n)
    {
        if (s.empty()) return false;

        Unsigned pos = 0;
        for (Unsigned i=0; i<n; i++)
        {
            pos = s.find("%s", pos);
            if (pos == std::string::npos) return false;
            pos += 2;
        }

        pos = s.find("%s", pos);
        if (pos != std::string::npos) return false;
            // only n of "%s " is allowed

        return true;
    }

    inline int cmd_pdfviewer(const std::string &s)
    {
        extern std::string PDF_VIEW;
        std::string msg;

        bool flag = test_format(s, 1);
        if (!flag)
            return cmd_err("Error : invalid format. The correct format is :\n"
            "\t viewername [params..]%s.pdf \n"
            "\n"
            "Examples :\n"
            "\t evince %s.pdf [default]\n"
            "\t xdg-open %s.pdf\n"
            "\t firefox %s.pdf [bad choice] \n");


        PDF_VIEW = s;
        msg += "[PDF_VIEW]        " + PDF_VIEW;
        return cmd_msg(msg);
    }

    inline int cmd_svg2pdf(const std::string &s)
    {
        extern std::string SVG_TO_PDF;

        bool flag = test_format(s, 2);
        if (!flag)
            return cmd_err("\n"
            "Error : invalid format. The correct format is :\n"
            "\t toolname [params]%s.pdf [params]%s.svg\n"
            "\n"
            "Examples :\n"
            "\t inkscape --export-pdf=%s.pdf %s.svg [default]\n"
            "\t rsvg-convert -f pdf -o %s.pdf %s.svg\n"
            );

        std::string msg;
        SVG_TO_PDF = s;
        msg += "[SVG_TO_PDF]        " + SVG_TO_PDF;
        return cmd_msg(msg);
    }

    inline int cmd_locale(const std::string &s)
    {
        std::string msg;
        LOCALE = s;
        msg += "[LOCALE]           " + LOCALE;
        cmd_msg(msg);
        return REINIT_NCURCES;
    }

    /*!
     * \brief change directory
     * \param s command parameters
     *      "" - list of used directories
     *      "WORK_DIR[W|w] dirname" - change work directory
     *      "OUTPUT_DIR[O|o] dirname" - change output directory
     */
    inline int cmd_change_dir(const std::string &s)
    {
        std::string filename = TMP_DIR + "/terminal.shot";
        std::string msg;

        auto test_and_assign = [](std::string s) -> std::string
        {
            auto pos = s.find_first_of(' ');
            std::string v = (pos == std::string::npos) ? "" : std::string(s, pos + 1);
            v = ltrim(v);

            if (v[0] == '~') return HOME_DIR + std::string(v, 1); // replace '~'
            if (std::strncmp(v.c_str(), "/", 1) == 0) return v;
            return HOME_DIR + '/' + v;
        };

        if (std::toupper(s[0]) == 'O')
        {
            OUTPUT_DIR = test_and_assign(s);
            return cmd_msg(OUTPUT_DIR);
        }

        if (std::toupper(s[0]) == 'W')
        {
            WORK_DIR = test_and_assign(s);
            errno = 0; int r = chdir(&WORK_DIR[0]);
            if (r != 0) return cmd_msg(strerror(errno));
            return cmd_msg(WORK_DIR);
        }
        return 0;
    }

    inline int cmd_env(const std::string &ss)
    {
        std::string s(ss);
        s = trim(s);

        if (s.empty())
        {
            std::string msg;
            msg += "\n";
            msg += "[LOCALE]           " + LOCALE + '\n';
            msg += "[WORK_DIR]         " + WORK_DIR + '\n';
            msg += "[OUTPUT_DIR]       " + OUTPUT_DIR + '\n';
            msg += "[PDF_VIEW]         " + PDF_VIEW + '\n';
            msg += "[SVG_TO_PDF]       " + SVG_TO_PDF + '\n';
            return cmd_msg(msg);
        }
        else
        {
            Unsigned pos = s.find_first_of(' ');
            std::string val = (pos == std::string::npos) ? "" : std::string(s, pos + 1);

            switch (std::toupper(s[0]))
            {
            case 'L' : return cmd_locale(val);
            case 'P' : return cmd_pdfviewer(val);
            case 'S' : return cmd_svg2pdf(val);
            case 'W' :
            case 'O' : return cmd_change_dir(s);
            }
        }
        return 0;
    }

    inline int cmd_cd(const std::string &filename)
    {
        errno = 0;
        int r = chdir(&filename[0]);
        std::string msg;

        if (r == 0)
        {
            char buffer[1024];
            char *wd = getcwd(buffer, 1024);
            WORK_DIR = wd;
            msg += wd;
        }
        else msg = strerror(errno);

        return cmd_msg(msg);
    }

    inline int cmd_flist(const std::string &)
    {
        std::string msg("\nfiles: \n\n");
        msg += "  st fileID       filename   rS   rD [list..]\n";
        msg += "---------------------------------------------\n";

        for (auto e : flowerdomes)
        {
            std::string status = " ";
            status += (read_flowerdome == e.first) ? "R" : "-";
            status += (write_flowerdome == e.first) ? "W" : "-";
            status += (iflowerdomes.count(e.first)) ? "^" : "-";

            std::string species = e.second->species_name;;

            Unsigned n = e.second->num_sprouts();
            Unsigned n1 = e.second->num_dynamic_sprouts();

            char buffer[80];
            std::sprintf(buffer, "%3s %6s %14s %3ld* %3ld  [%s]",
                status.c_str(), e.first.c_str(), species.c_str(), n-n1, n1,
                e.second->get_info().c_str()
            );
            msg += std::string(buffer) + "\n";
        }
        msg += "\n";
        if (read_flowerdome == "") msg += "no file with R status\n";
        if (write_flowerdome == "") msg += "no file with W status\n";

        return cmd_msg(msg);
    }

    inline int cmd_fadd(const std::string &s)
    {
        if (s.empty()) return cmd_msg("\nError: no fileID");

        std::vector<std::string> ss = split(s, ' ');

        std::string name = ss[0];
        flowerdomes[name] = shared<Flowerdome>(new Flowerdome);

        std::string status = (ss.size() > 1) ? ss[1] : "";
        for (auto &ch : status) ch = std::toupper(ch);

        write_flowerdome = name;
        wflowerdome = flowerdomes[name];

        return cmd_flist("");
    }

    inline int cmd_fdel(const std::string &fileID)
    {
        if (fileID.empty()) return cmd_msg("\nError: no fileID");

        if (read_flowerdome == fileID)
        {
            read_flowerdome = "";
            rflowerdome = nullptr;
        }
        if (write_flowerdome == fileID)
        {
            write_flowerdome = "";
            wflowerdome = nullptr;
        }
        flowerdomes.erase(fileID);
        iflowerdomes.erase(fileID);

        return cmd_flist("");
    }

    inline int cmd_fstatus(const std::string &param)
    {
        if (param.empty()) return cmd_flist("");

        auto tokens = split(param, ' ');

        std::string fileID = tokens[0];
        std::string status = (tokens.size() > 1) ? tokens[1] : "";
        for (auto &ch : status) ch = std::toupper(ch);

        if (status == "W")
        {
            write_flowerdome = fileID;
            wflowerdome = flowerdomes[fileID];
        }
        else if (status == "R")
        {
            read_flowerdome = fileID;
            rflowerdome = flowerdomes[fileID];
        }
        else if (status == "RW" || status == "WR")
        {
            read_flowerdome = write_flowerdome = fileID;
            wflowerdome = flowerdomes[fileID];
            rflowerdome = flowerdomes[fileID];
        }
        else status = "";

        return cmd_flist("");
    }

    inline int cmd_fload(const std::string &s)
    {
        if (s.empty()) return cmd_err("\nError: no filename");

        std::vector<std::string> ss = split(s, ' ');

        std::string filename = ss[0];
        std::string ID = (ss.size() > 1) ? ss[1] : "noname";

        flowerdomes[ID] = shared<Flowerdome>(new Flowerdome);
        std::string msg = flowerdomes[ID]->load(filename);

        read_flowerdome = ID;
        rflowerdome = flowerdomes[ID];

        if (write_flowerdome.empty())
        {
            write_flowerdome = ID;
            wflowerdome = flowerdomes[ID];
        }

        return cmd_flist("");
    }

    inline int cmd_fsave(const std::string &s)
    {
        if (s.empty()) return cmd_err("\nError: no params");

        if (dir_exists(OUTPUT_DIR.c_str()) != 1)
        {
            std::string cmd("mkdir ");
            cmd += OUTPUT_DIR;
            SYSTEM(cmd);
        }

        std::string filename;
        std::vector<std::string> tokens = split(s, ' ');
        if (tokens.size() == 1) return cmd_err("\nError: false format");

        std::string ID = tokens[0];
        filename = OUTPUT_DIR + "/" + tokens[1];
        flowerdomes[ID]->save(filename);

        return cmd_msg(filename + " saved");
    }

    inline int cmd_fimprove(const std::string &param)
    {
        if (param.empty()) return cmd_msg("Error: no params");

        auto tokens = split(param);
        std::string method = tokens[0];
        std::string fileID = (tokens.size() >=2) ? tokens[1] : read_flowerdome;

        if (method == "MLM")
        {
            flowerdomes[fileID]->recover_missed_data_by_direct_MLM();
            iflowerdomes.insert(fileID);
            return cmd_flist("");
        }
        else if (method == "random")
        {
            flowerdomes[fileID]->recover_missed_data_by_random();
            iflowerdomes.erase(fileID);
            return cmd_flist("");
        }

        //if (std::tolower(param[0]) == 'b') rflowerdome->recover_missed_data_by_direct_MLM(BUD);
        //if (std::tolower(param[0]) == 'f') rflowerdome->recover_missed_data_by_direct_MLM(FLOWER);
        return cmd_msg("Error: incorrect params");
    }

    inline int cmd_mlist(const std::string &)
    {
        std::string msg("\nmodels : \n\n");
        msg += " st modelID fileID:sproutID\n";
        msg += "---------------------------\n";
        for (auto e : sproutmodels)
        {
            std::string status = " ";
            if (e.first == active_model) status = " R";

            char buffer[80];
            std::sprintf(buffer, " %2s %7s %15s ",
                status.c_str(), e.first.c_str(),
                e.second->get_info().c_str()
            );
            msg += std::string(buffer) + "\n";
        }
        return cmd_msg(msg);
    }

    inline int cmd_mstatus(const std::string &params)
    {
        if (params.empty())
            return cmd_mlist("");

        std::vector<std::string> tokens = split(params);
        std::string modelID = tokens[0];

        if (sproutmodels.count(modelID) == 0)
            return cmd_err("\nError : no model with ID given.");

        active_model = modelID;
        sproutmodel = sproutmodels[active_model];

        return cmd_mlist("");
    }

    inline int cmd_madd(const std::string &params)
    {
        if (params.empty())
            return cmd_err("\nError: no fileID");

        std::vector<std::string> tokens = split(params);
        std::string fileID_and_sproutID = tokens[0];
        std::string modelID = tokens.size() > 1 ? tokens[1] : "noname";

        std::vector<std::string> ids = split(fileID_and_sproutID, ':');
        if (ids.size() != 2) return cmd_err("Error: false fileID:sproutID format");

        std::string fileID = ids[0];
        Unsigned sproutID = std::atol(ids[1].c_str());

        if (flowerdomes.count(fileID) == 0)
            return cmd_err("Error : no file with ID given.");

        if (!flowerdomes[fileID]->has_sprout(sproutID))
            return cmd_err("Error : no sprout with ID given.");

        std::string info = fileID + ":" + ids[1];
        sproutmodel = shared<Sproutmodel>(new Sproutmodel);
        sproutmodel->build(*flowerdomes[fileID], sproutID, info);
        sproutmodels[modelID] = sproutmodel;
        active_model = modelID;

        return cmd_mlist("");
    }

    inline int cmd_mdel(const std::string &params)
    {
        if (sproutmodels.empty())
            return cmd_err("\nError : no models");

        if (params.empty())
            return cmd_err("\nError : no modelID.");

        std::vector<std::string> tokens = split(params);
        std::string modelID = tokens[0];

        if (sproutmodels.count(modelID) == 0)
            return cmd_err("\nError : no model with modelID given.");

        sproutmodels.erase(modelID);
        if (active_model == modelID) active_model = "";
        if (!sproutmodels.empty())
        {
            active_model = sproutmodels.begin()->first;
            sproutmodel = sproutmodels[active_model];
        }
        else sproutmodel = nullptr;

        return cmd_mlist("");
    }

    inline int cmd_plots(const std::string &params)
    {
        std::string msg;

        if (sproutmodel == nullptr)
            return cmd_err("\nError: no models");

        if (params.empty())
            return cmd_msg(sproutmodel->view_list());

        std::string s(params);
        s = trim(s);

        Unsigned id = std::atoi(s.c_str());
        active_plot = id;

        if (id >= sproutmodel->view_list_size())
            return cmd_err("\nError : incorrect plotID\n");

        sproutmodel->view(id);
        return cmd_msg("\n Done");
    }

    inline int cmd_padd(const std::string &params)
    {
        std::string msg;

        if (sproutmodel == nullptr)
            return cmd_err("\nError: no models yet");

        if (params.empty()) return cmd_msg(sproutmodel->view_list());

        double width = default_width;
        auto pos = params.find_first_of('w');
        if (pos != std::string::npos)
        {
            std::string s(params, pos + 1);
            s = trim(s);
            width = std::atof(s.c_str());
        }

        std::string s(params, 0, pos);
        s = trim(s);
        if (!s.empty()) active_plot = std::atoi(s.c_str());

        dependence_item &d = sproutmodel->dependencies_items.at(active_plot);

        dependence_item::PF2 pf2 = d.pf2;
        ((*sproutmodel).*pf2)(d, width);

        gnuplotParams v = d.gnuparam;
        plot(v);

        return cmd_msg("\n Done");
    }

    inline int cmd_pdel(const std::string &params)
    {
        std::string msg;

        if (sproutmodel == nullptr)
            return cmd_err("\nError: no models yet");

        if (!params.empty())
        {
            std::string s(params);
            s = trim(s);

            Unsigned id = std::atoi(s.c_str());
            Unsigned n = sproutmodel->view_list_size();
            if (id >= n) return cmd_err("\nError : incorrect plotID");
            active_plot = id;
        }

        dependence_item &d = sproutmodel->dependencies_items.at(active_plot);

        gnuplotParams &v = d.gnuparam;
        v.clear();
        plot(v);

        return cmd_msg("\n Done");
    }

    inline int cmd_test(const std::string &params)
    {
        std::string msg;

        if (sproutmodel == nullptr)
            return cmd_err("\nError: no models");

        if (params.empty()) return cmd_err("\nError: no params");

        auto tokens = split(params);

        Unsigned id = std::atol(tokens[0].c_str());
        if (id >= sproutmodel->dependencies_items.size()) return cmd_err("Error: incoorect plotID");
        dependence_item &d = sproutmodel->dependencies_items.at(id);
        Unsigned data_n = d.normality_test_n;

        Unsigned n = tokens.size() == 1 ? data_n : std::atol(tokens[1].c_str());

        dependence_item::PF3 pf3 = d.pf3;
        msg = ((*sproutmodel).*pf3)(d, n);

        return cmd_msg(msg);
    }

    inline int cmd_scheme(const std::string &param)
    {
        if (param.empty()) return cmd_err("\nError: no params");
        if (sproutmodels.empty()) return cmd_err("\nError: no models");

        auto tokens = split(param);
        std::string date = tokens.size() == 2 ? tokens[1] : "0";

        std::string modelID = "";
        std::string fileID = "";
        Unsigned sproutID = 0;

        auto toks = split(tokens[0], ':');
        if (toks.size() == 1)
        {
            modelID = toks[0];
            if (sproutmodels.count(modelID) == 0)
                return cmd_err("\nError: incorrect <modelID>");

            std::string fileID_and_sproutID = sproutmodels[modelID]->get_info();

            std::vector<std::string> ids = split(fileID_and_sproutID, ':');
            if (ids.size() != 2) return cmd_err("Error: false fileID:sproutID format");

            fileID = ids[0];
            sproutID = std::atol(ids[1].c_str());
        }
        else
        {
            fileID = toks[0];
            sproutID = std::atol(toks[1].c_str());
        }

        if (flowerdomes.count(fileID) == 0) return cmd_err("\nError: no file <fileID>");

        std::string msg;
        std::string svgname = TMP_DIR + "/terminal";
        std::string pdfname = TMP_DIR + "/terminal";

        Sprout &sprout = flowerdomes[fileID]->sprout(sproutID);

        Time time;
        if (is_date(date)) time = Time(date.c_str());
        else
        {
            time = sprout.day_of_first_start(FLOWER);
            time += atoi(date.c_str());
        }

        SproutScheme scheme = sprout.scheme_generate(time);
        sprout.scheme_print(svgname + ".svg", scheme);

        char buffer[1024];
        std::sprintf(buffer, SVG_TO_PDF.c_str(), pdfname.c_str(), svgname.c_str());
        SYSTEM(buffer);

        return cmd_msg("\n Done");
    }

    inline int cmd_video(const std::string &param)
    {
        if (param.empty()) return cmd_err("\nError: no params");
        if (sproutmodels.empty()) return cmd_err("\nError: no models");

        auto tokens = split(param);
        if (tokens.size() != 2) return cmd_err("\nError: incorrect format");

        std::string modelID = tokens[0];
        if (sproutmodels.count(modelID) == 0) return cmd_err("\nError: incorrect <modelID>");

        std::string fileID_and_sproutID = sproutmodels[modelID]->get_info();

        std::vector<std::string> ids = split(fileID_and_sproutID, ':');
        if (ids.size() != 2) return cmd_err("Error: false fileID:sproutID format");

        std::string fileID = ids[0];
        Unsigned sproutID = std::atol(ids[1].c_str());

        if (flowerdomes.count(fileID) == 0) return cmd_err("\nError: no file <fileID>");

        std::string dates = tokens[1];
        ids = split(dates, ':');
        if (ids.size() != 2) return cmd_err("Error: false <date1:date2> format");

        Sprout &sprout = flowerdomes[fileID]->sprout(sproutID);

        std::string msg;
        std::string filename = TMP_DIR + "/terminal";

        Time time1, time2;
        if (is_date(ids[0]))
        {
            time1 = Time(ids[0].c_str());
            time2 = Time(ids[1].c_str());
        }
        else
        {
            time1 = sprout.day_of_first_start(FLOWER);
            time1 += atoi(ids[0].c_str());

            time2 = sprout.day_of_first_start(FLOWER);
            time2 += atoi(ids[1].c_str());
        }

        for (Time t=time1; t<=time2; ++t)
        {
            SproutScheme scheme = sprout.scheme_generate(t);
            std::string caption = t.to_string();
            sprout.scheme_print(filename + ".svg", scheme, caption);

            char buffer[1024];
            std::sprintf(buffer, SVG_TO_PDF.c_str(), filename.c_str(), filename.c_str());
            SYSTEM(buffer);

            std::this_thread::sleep_for(1s);
        }

        return cmd_msg("\n Done");
    }

    /*!
     * \brief open graphical terminal to view current result
     */
    inline int cmd_window(const std::string &)
    {
        std::string filename = TMP_DIR + "/terminal";

        struct stat sb;
        if (stat((filename + ".pdf").c_str(), &sb) == -1)
        {
            Dimensions dimensions(400, 600);
            Document doc(filename + ".svg", Layout(dimensions, Layout::TopLeft));
            doc << Text(Point(40,20), "Wellcome to ornamental", Fill(Color::Blue),
                 Font(12, "Verdana"), Stroke(0.1, Color::Blue));
            doc.save();

            char buffer[1024];
            std::sprintf(buffer, SVG_TO_PDF.c_str(), filename.c_str(), filename.c_str());
            SYSTEM(buffer);
        }

        if (PDF_VIEW.empty())
            return cmd_err("\nError : set PDF_VIEW");

        auto pid = fork();
        if (pid == 0)
        {
            char buffer[1024];
            std::sprintf(buffer, PDF_VIEW.c_str(), filename.c_str());
            SYSTEM(buffer);
            exit(0);
        }
        return cmd_msg("\n Done");
    }

    inline int cmd_psave(const std::string &s)
    {
        if (s.empty()) return cmd_err("\nError: no params");

        if (dir_exists(OUTPUT_DIR.c_str()) != 1)
        {
            std::string cmd("mkdir ");
            cmd += OUTPUT_DIR;
            SYSTEM(cmd);
        }

        std::string filename;
        filename = OUTPUT_DIR + "/" + s;
        std::string cmd("cp ");
        cmd += std::string(getenv("HOME")) + "/.ornamental/terminal.pdf ";
        cmd += OUTPUT_DIR + "/" + s;
        SYSTEM(cmd);

        return cmd_msg(filename + " saved");
    }

    inline int cmd_gen(const std::string &params)
    {
        // generate|g <modelID> N <n1=0:n2=*>

        if (params.empty()) return cmd_err("\nError: no params");
        auto tokens = split(params);
        if (tokens.size() < 2) return cmd_err("\nError: not enough params");

        std::string modelID = tokens[0];
        if (sproutmodels.count(modelID) == 0) return cmd_err("\nError: no model");

        Unsigned n = std::atol(tokens[1].c_str());
        std::string n1n2 = (tokens.size() == 3)? tokens[2] : "*:*";

        auto toks = split(n1n2, ':');
        if (toks.size() == 1) return cmd_err("\nError: incorrect format");
        Unsigned n1 = (toks[0] == "*") ? inf<Unsigned>() : std::atol(toks[0].c_str());
        Unsigned n2 = (toks[1] == "*") ? inf<Unsigned>() : std::atol(toks[1].c_str());

        auto model = sproutmodels[modelID];

        Unsigned last_id = wflowerdome->empty() ? 0 : wflowerdome->rbegin()->first;
        for (Unsigned i=1; i<=n; i++)
        {
            shared<Sprout> sprout = model->generate_sprout(n1, n2);
            (*wflowerdome)[last_id + i] = *sprout;
        }

        std::string msg = "\n" + std::to_string(n) + " model sprouts added to \'"
            + write_flowerdome + "\'\n";
        return cmd_msg(msg);
    }

#define ITEM(s1, s2, cmd)  {s1, &cmd}, {s2, &cmd},
    const std::map<std::string, PF> orn_function_table =
    {
    ITEM("help", "h", cmd_help)
    ITEM("environment", "env", cmd_env)
    ITEM("window", "w", cmd_window)
    ITEM("chdir", "cd", cmd_cd)

    ITEM("files", "f", cmd_fstatus)
    ITEM("fload", "f<", cmd_fload)
    ITEM("fsave", "f>", cmd_fsave)
    ITEM("file^", "f^", cmd_fimprove)
    ITEM("fadd", "f+", cmd_fadd)
    ITEM("fdel", "f-", cmd_fdel)

    ITEM("models", "m", cmd_mstatus)
    ITEM("madd", "m+", cmd_madd)
    ITEM("mdel", "m-", cmd_mdel)

    ITEM("plots", "p", cmd_plots)
    ITEM("padd", "p+", cmd_padd)
    ITEM("pdel", "p-", cmd_pdel)
    ITEM("psave", "p>", cmd_psave)

    ITEM("scheme", "scheme", cmd_scheme)
    ITEM("video", "video", cmd_video)
    ITEM("test", "t", cmd_test)
    ITEM("generate", "g", cmd_gen)
    };
#undef ITEM

}
#endif // ORN_COMMAND_H
