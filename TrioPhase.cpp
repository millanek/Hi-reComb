//
//  TrioPhase.cpp
//  Hi-reComb
//
//  Created by Milan Malinsky on 23.08.23.
//

#include "TrioPhase.hpp"

#define SUBPROGRAM "TrioPhase"

#define DEBUG 1

static const char *TRIO_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] OFFSPRING.vcf(.gz) PARENTS.vcf(.gz) PARENT1,PARENT2\n"
"Calculate the Allele Frequencies per population/species from a VCF \n"
"\n"
HelpOption RunNameOption
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "hn:";


static const struct option longopts[] = {
    { "help",   no_argument, NULL, 'h' },
    { "run-name",   required_argument, NULL, 'n' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string offspringVcfFile;
    static string parentsVcfFile;
    static string runName = "out";
    static string parentNames;
}


int trioPhaseMain(int argc, char** argv) {
    parseTrioPhaseOptions(argc, argv);
    std::vector<string> parentNames = split(opt::parentNames, ',');
    if (parentNames.size() != 2) {
        std::cerr << "ERROR: The parents need to be exactly two individual IDs separated by a comma. Exiting..." << std::endl; exit(1);
    }
        
    std::vector<std::string> fields;
    string chr; string coord;
    std::map<int,string> offspringPosToGT;
    // std::ostream* outFileAF = createWriter(stripExtension(opt::setsFile) + "_" + opt::runName + "_AF" + ".txt");
    
    string line; // for reading the input files
    std::istream* offspringVcfFile = createReader(opt::offspringVcfFile.c_str());
    std::istream* parentsVcfFile = createReader(opt::parentsVcfFile.c_str());
    std::ostream* outFilePhasedHets = createWriter(opt::runName + "_phasedHets.txt");
    *outFilePhasedHets << "BLOCK:" << std::endl;
    std::cerr << "INFO: Loading the offspring VCF file...";
    while (getline(*offspringVcfFile, line)) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end()); // Deal with any left over \r from files prepared on Windows
        if (line[0] == '#' && line[1] == '#')
            continue;
        else if (line[0] == '#' && line[1] == 'C') {
            fields = split(line, '\t');
            std::vector<std::string> sampleNames(fields.begin()+NUM_VCF_NON_GENOTYPE_COLUMNS,fields.end());
            if (sampleNames.size() > 1) {
                std::cerr << "ERROR: The offspring VCF file should only contain one sample." << std::endl; exit(1);
            }
        } else {
            fields = split(line, '\t');
            VariantInfo v(fields); if (v.onlyIndel) continue; // Only consider biallelic SNPs
            if (v.refAllele.length() > 1) continue;
            if (v.altAlleles.size() > 1) continue;
            if (v.altAlleles[0].length() > 1) continue;
            
            std::vector<std::string> genotypes(fields.begin()+NUM_VCF_NON_GENOTYPE_COLUMNS,fields.end());
            
            
            string firstAllele; firstAllele += getAllele(genotypes[0], 0, v);
            string secondAllele; secondAllele += getAllele(genotypes[0], 2, v);
            
            if (firstAllele != secondAllele) { // Only recording hets
                offspringPosToGT[v.posInt] = firstAllele + secondAllele;
            }
            
            genotypes.clear(); genotypes.shrink_to_fit();
        }
    }
    std::cerr << "DONE" << std::endl;
    std::cerr << "INFO: Number of hets: " << offspringPosToGT.size() << std::endl;
    
    std::clock_t startTime = std::clock(); int totalVariantNumber = 0;
    
    int numViolations = 0; int numPhasedHets = 0; int numParentsUninformative = 0; int remaining = 0;
    int iP1; int iP2;
    std::cerr << "INFO: Processing the parents VCF file...";
    while (getline(*parentsVcfFile, line)) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end()); // Deal with any left over \r from files prepared on Windows
        if (line[0] == '#' && line[1] == '#')
            continue;
        else if (line[0] == '#' && line[1] == 'C') {
            fields = split(line, '\t');
            std::vector<std::string> sampleNames(fields.begin()+NUM_VCF_NON_GENOTYPE_COLUMNS,fields.end());
            
            std::vector<std::string>::iterator itP1 = std::find(sampleNames.begin(), sampleNames.end(), parentNames[0]);
            std::vector<std::string>::iterator itP2 = std::find(sampleNames.begin(), sampleNames.end(), parentNames[1]);
            if (itP1 == sampleNames.end()) {
                std::cerr << "ERROR: Could not locate parent 1 in the VCF. Exiting..." << std::endl; exit(1);
            }
            if (itP2 == sampleNames.end()) {
                std::cerr << "ERROR: Could not locate parent 2 in the VCF. Exiting..." << std::endl; exit(1);
            }
            
            iP1 = itP1 - sampleNames.begin(); iP2 = itP2 - sampleNames.begin();
        } else {
            totalVariantNumber++;
           // if (totalVariantNumber % reportProgressEvery == 0) reportProgessVCF(totalVariantNumber, startTime);
            
            fields = split(line, '\t');
            VariantInfo v(fields); if (v.onlyIndel) continue; // Only consider SNPs
            
            string offspringGT;
            if (offspringPosToGT.count(v.posInt) == 1) {
                offspringGT = offspringPosToGT.at(v.posInt);
                char offspringAllele1 = offspringGT[0];
                char offspringAllele2 = offspringGT[1];
                
                std::vector<std::string> genotypes(fields.begin()+NUM_VCF_NON_GENOTYPE_COLUMNS,fields.end());
                if (genotypes[iP1][0] == '.' || genotypes[iP1][2] == '.') continue;
                if (genotypes[iP2][0] == '.' || genotypes[iP2][2] == '.') continue;
                
                char P1firstAllele = getAllele(genotypes[iP1], 0, v);
                char P1secondAllele = getAllele(genotypes[iP1], 2, v);
                
                // Offspring: GC P1: AA
                // If neither offspring genotype matches Parent 1, we classify this SNP as a Mendelian violation and cannot phase:
                if (offspringAllele1 != P1firstAllele && offspringAllele1 != P1secondAllele && offspringAllele2 != P1firstAllele && offspringAllele2 != P1secondAllele) { numViolations++; continue; }
                
                char P2firstAllele = getAllele(genotypes[iP2], 0, v);
                char P2secondAllele = getAllele(genotypes[iP2], 2, v);
                
                // If neither offspring genotype matches Parent 1, we classify this SNP as a Mendelian violation and cannot phase:
                if (offspringAllele1 != P2firstAllele && offspringAllele1 != P2secondAllele && offspringAllele2 != P2firstAllele && offspringAllele2 != P2secondAllele) { numViolations++; continue; }
                
                // If offspringAllele1 is not present in either parent, we classify this SNP as a Mendelian violation and cannot phase:
                if (offspringAllele1 != P1firstAllele && offspringAllele1 != P1secondAllele && offspringAllele1 != P2firstAllele && offspringAllele1 != P2secondAllele) { numViolations++; continue; }
                
                // If offspringAllele2 is not present in either parent, we classify this SNP as a Mendelian violation and cannot phase:
                if (offspringAllele2 != P1firstAllele && offspringAllele2 != P1secondAllele && offspringAllele2 != P2firstAllele && offspringAllele2 != P2secondAllele) { numViolations++; continue; }
                
                // If both parents are heterozygous, we cannot phase
                // Offspring: 0/1
                // Parents P1: 0/1 P2: 0/1    -->   Impossible to phase
                if (P1firstAllele != P1secondAllele && P2firstAllele != P2secondAllele) { numParentsUninformative++; continue; }
                
                // Offspring: 0/1
                // Parents P1: 0/0 P2: 1/1
                // Parents P1: 1/1 P2: 0/0
                if (P1firstAllele == P1secondAllele && P2firstAllele == P2secondAllele) {
                    if (offspringAllele1 == P1firstAllele) { // Parents P1: 0/0 P2: 1/1
                        numPhasedHets++; printPhasedLine(*outFilePhasedHets, numPhasedHets, v, offspringAllele1, offspringAllele2); continue;
                    }
                    if (offspringAllele1 == P2firstAllele) { // Parents P1: 1/1 P2: 0/0
                        numPhasedHets++; printPhasedLine(*outFilePhasedHets, numPhasedHets, v, offspringAllele2, offspringAllele1);
                        continue;
                    }
                }
                
                // Offspring: A/T
                // Parents P1:AT P2:TT
                // Parents P1:AT P2:AA
                if (P1firstAllele != P1secondAllele && P2firstAllele == P2secondAllele) {
                    if (offspringAllele2 == P2firstAllele) { // Parents P1:AT P2:TT
                        numPhasedHets++; printPhasedLine(*outFilePhasedHets, numPhasedHets, v, offspringAllele1, offspringAllele2); continue;
                    }
                    if (offspringAllele1 == P2firstAllele) { // Parents P1:TT P2:AT
                        numPhasedHets++; printPhasedLine(*outFilePhasedHets, numPhasedHets, v, offspringAllele2, offspringAllele1); continue;
                    }
                }
                
                // Offspring: A/T
                // Parents P1:AA P2:AT
                // Parents P1:TT P2:AT
                if (P1firstAllele == P1secondAllele && P2firstAllele != P2secondAllele) {
                    if (offspringAllele1 == P1firstAllele) { // Parents P1:AA P2:AT
                        numPhasedHets++; printPhasedLine(*outFilePhasedHets, numPhasedHets, v, offspringAllele1, offspringAllele2); continue;
                    }
                    if (offspringAllele2 == P1firstAllele) { // Parents P1:TT P2:AT
                        numPhasedHets++; printPhasedLine(*outFilePhasedHets, numPhasedHets, v, offspringAllele2, offspringAllele1); continue;
                    }
                }
                
                std::cout << "offspringAllele1: " << offspringAllele1 << "\toffspringAllele2: " << offspringAllele2 << std::endl;
                std::cout << "P1firstAllele: " << P1firstAllele << "\tP1secondAllele: " << P1secondAllele << std::endl;
                std::cout << "P2firstAllele: " << P2firstAllele << "\tP2secondAllele: " << P2secondAllele << std::endl;
                
                genotypes.clear(); genotypes.shrink_to_fit();
            }
        }
    }
    std::cerr << "DONE" << std::endl;
    
    int parentsNotPolymorphic = offspringPosToGT.size() - numViolations - numPhasedHets - numParentsUninformative - remaining;
    std::cout << "INFO: Number of phased hets: " << numPhasedHets << " (" << (double)numPhasedHets/offspringPosToGT.size() << ")" << std::endl;
    std::cout << "INFO: Number of unphased hets (parents uninformative): " << numParentsUninformative << " (" << (double)numParentsUninformative/offspringPosToGT.size() << ")" << std::endl;
    std::cout << "INFO: Not polymorphic in parents: " << parentsNotPolymorphic << " (" << (double)parentsNotPolymorphic/offspringPosToGT.size() << ")" << std::endl;
    std::cout << "INFO: Number of mendelian violations: " << numViolations << " (" << (double)numViolations/offspringPosToGT.size() << ")" << std::endl;
    if (remaining > 0) std::cout << "INFO: Unknown: " << remaining << std::endl;
    
    return 0;
}



void parseTrioPhaseOptions(int argc, char** argv) {
    bool die = false; string regionArgString; std::vector<string> regionArgs;
    std::vector<string> windowSizeStep;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'n': arg >> opt::runName; break;
            case 'h':
                std::cout << TRIO_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (argc - optind < 3) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > 3)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << TRIO_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::offspringVcfFile = argv[optind++];
    opt::parentsVcfFile = argv[optind++];
    opt::parentNames = argv[optind++];
}
