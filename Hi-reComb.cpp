//
//  main.cpp
//  Hi-reComb
//
//  Created by Milan Malinsky on 02.09.22.
//

// Hi-reComb

// Some questions:
// 1) Are insert length distributions along the genome random in the joint recomb + non-recomb pairs?
// 2) Does the average genome-wide recombination rate inferred depend on the insert sizes that are used? E.g., try inferring with only 0-1000bp, 1000-2000bp, 2000-10000bp, and longer? Specifically, at less than 1000bp, maybe there will be an effect of gene-conversion?
// 3) Find the coverage of informative read-pairs at every covered het SNP.

#include "generalUtils.hpp"
#include "recombFromInformativePairsSAM.hpp"
#include "findInformativePairs.hpp"
#include "countMendelianViolations.hpp"
#include "TrioPhase.hpp"


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
"Extra tools:\n"
"           CountViolations     \n"
"           TrioPhase     \n"
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
