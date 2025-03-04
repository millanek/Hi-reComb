//
//  main.cpp
//  Hi-reComb
//
//  Created by Milan Malinsky on 02.09.22.
//

// Hi-reComb

// Some questions:
// 1) Are insert length distributions along the genome random in the joint recomb + non-recomb pairs?

#include "generalUtils.hpp"
#include "recombFromInformativePairsSAM.hpp"
#include "findInformativePairs.hpp"
#include "countMendelianViolations.hpp"
#include "TrioPhase.hpp"
#include "simulateAndReconstruct.hpp"


#define AUTHOR "Milan Malinsky"
#define PACKAGE_VERSION "0.1 r1"


static const char *VERSION_MESSAGE =
"evo software Version " PACKAGE_VERSION "\n"
"Written by Milan Malinsky.\n"
"\n";

static const char *USAGE_MESSAGE =
"Program: " PROGRAM_BIN "\n"
"Version: " PACKAGE_VERSION "\n"
"Contact: " AUTHOR " [" PACKAGE_BUGREPORT "]\n"
"Usage: " PROGRAM_BIN " <command> [options]\n\n"
"Commands:\n"
"           FindInfoPairs       Find read pairs that would be informative for estimating the recombination\n"
"           RecombMap           Estimate the recombination map from informative Hi-C read pairs\n"
"Utilities:\n"
"           Simulate            Simulate informative read-pairs for a given map to evaluate confidence in reconstruction\n"
"           TrioPhase           Generate a phased het file for use with Hi-Recomb from a VCF with trio(s) (mother-father-offspring)\n"
// "           CountViolations     (specific use case) Find trios in a VCF\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

int main(int argc, char **argv) {
    
    if(argc <= 1)
    {
        std::cout << USAGE_MESSAGE;
        return 0;
    }
    else
    {
        std::string command(argv[1]);
        if(command == "help" || command == "--help" || command == "-h")
        {
            std::cout << USAGE_MESSAGE;
            return 0;
        }
        else if(command == "version" || command == "--version")
        {
            std::cout << VERSION_MESSAGE;
            return 0;
        }
        
        if(command == "RecombMap")
            RecombFromSAMMain(argc - 1, argv + 1);
        else if(command == "FindInfoPairs")
            InfoReadsMain(argc - 1, argv + 1);
        else if (command == "Simulate")
            SimulationMain(argc - 1, argv + 1);
        else if (command == "CountViolations")
            vioMain(argc - 1, argv + 1);
        else if (command == "TrioPhase")
            trioPhaseMain(argc - 1, argv + 1);
        else
        {
            std::cerr << "Unrecognized command: " << command << "\n";
            return 1;
        }
        return 0;
    }
}
