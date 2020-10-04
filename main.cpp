#include "orn/command.h"

namespace prg
{
    //---------------------------------------------------------------------------------------------
    //      global variables that are initialized manually or by an initialization file
    //---------------------------------------------------------------------------------------------
    std::string HOME_DIR = "";
    std::string TMP_DIR = "";
    std::string WORK_DIR = "";
    std::string OUTPUT_DIR = "";
    std::string LOCALE = "";
    std::string PDF_VIEW = "";
    std::string SVG_TO_PDF = "";

    shared<Flowerdome> rflowerdome; // current file to read
    shared<Flowerdome> wflowerdome; // current file to write
    shared<Sproutmodel> sproutmodel; // current model


    Unsigned active_sprout = 0; // current sprout to build model
    std::string active_model = ""; // current model to plot
    Unsigned active_plot = 0; // current plot

    std::map<std::string, shared<Flowerdome>> flowerdomes; // all files
    std::map<std::string, shared<Sproutmodel>> sproutmodels; // all models
    std::set<std::string> iflowerdomes; // list of improved flowerdomes

    std::string read_flowerdome = ""; // current file ID to read
    std::string write_flowerdome = ""; // current file ID to write
    std::string read_model = ""; // current model ID to read
    std::string write_model = ""; // current model ID to write
}

int main(int argc, char *argv[])
{
    set_ranseed(0); // debug mode

    const char *init_filename = (argc == 2) ? argv[1] : "";
    prg::cmd_Queue queue(&prg::orn_function_table, init_filename);
    queue.cycle();
}
