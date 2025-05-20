//
//  countMendelianViolations.cpp
//  Hi-reComb
//
//  Created by Milan Malinsky on 22.08.23.
//

#include "countMendelianViolations.hpp"

#define SUBPROGRAM "CountViolations"

#define DEBUG 1

static const char *VIO_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] OFFSPRING.vcf(.gz) PARENTS.vcf(.gz)\n"
"Calculate the Allele Frequencies per population/species from a VCF \n"
"\n"
HelpOption RunNameOption
"       -g, --use-genotype-probabilities        (optional) use genotype probabilities (GP tag) if present\n"
"                                               if GP not present calculate genotype probabilities from likelihoods (GL or PL tags) using a Hardy-Weinberg prior\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "hn:";


static const struct option longopts[] = {
    { "help",   no_argument, NULL, 'h' },
    { "run-name",   required_argument, NULL, 'n' },
    { "use-genotype-probabilities", no_argument, NULL, 'g'},
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string offspringVcfFile;
    static string parentsVcfFile;
    static string runName = "out";
    static bool useGenotypeProbabilities = false;
}


int vioMain(int argc, char** argv) {
    parseVioOptions(argc, argv);
        
    vector<string> fields;
    int reportProgressEvery = 10000; string chr; string coord;
    map<int,string> offspringPosToGT;
   // std::ostream* outFileAF = createWriter(stripExtension(opt::setsFile) + "_" + opt::runName + "_AF" + ".txt");
    
    string line; // for reading the input files
    std::istream* offspringVcfFile = createReader(opt::offspringVcfFile.c_str());
    std::istream* parentsVcfFile = createReader(opt::parentsVcfFile.c_str());
    std::cerr << "INFO: Loading the offspring VCF file" << std::endl;
    while (getline(*offspringVcfFile, line)) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end()); // Deal with any left over \r from files prepared on Windows
        if (line[0] == '#' && line[1] == '#')
            continue;
        else if (line[0] == '#' && line[1] == 'C') {
            fields = split(line, '\t');
            vector<string> sampleNames(fields.begin()+NUM_VCF_NON_GENOTYPE_COLUMNS,fields.end());
            if (sampleNames.size() > 1) {
                std::cerr << "ERROR: The offspring VCF file should only contain one sample." << std::endl; exit(1);
            }
        } else {
            fields = split(line, '\t');
            VariantInfo v(fields); if (v.onlyIndel) continue; // Only consider biallelic SNPs
            if (v.refAllele.length() > 1) continue;
            if (v.altAlleles.size() > 1) continue;
            if (v.altAlleles[0].length() > 1) continue;
            
            vector<string> genotypes(fields.begin()+NUM_VCF_NON_GENOTYPE_COLUMNS,fields.end());
            
            string firstAllele; string secondAllele;
            if (genotypes[0][0] == '0') firstAllele = v.refAllele[0];
            else firstAllele = v.altAlleles[0][0];
            if (genotypes[0][2] == '0') secondAllele = v.refAllele[0];
            else secondAllele = v.altAlleles[0][0];
            
            offspringPosToGT[v.posInt] = firstAllele + secondAllele;
 
            // Only consider biallelic SNPs
            genotypes.clear(); genotypes.shrink_to_fit();
        }
    }
    std::cerr << "INFO: DONE..." << std::endl;
    
    std::clock_t startTime = std::clock(); int totalVariantNumber = 0;
    
    vector<int> numViolations; int numSharedBetweenVCFs = 0;
    std::cerr << "INFO: Processing the parents VCF file.." << std::endl;
    while (getline(*parentsVcfFile, line)) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end()); // Deal with any left over \r from files prepared on Windows
        if (line[0] == '#' && line[1] == '#')
            continue;
        else if (line[0] == '#' && line[1] == 'C') {
            fields = split(line, '\t');
            vector<string> sampleNames(fields.begin()+NUM_VCF_NON_GENOTYPE_COLUMNS,fields.end());
            print_vector(sampleNames, std::cout);
            numViolations.resize(sampleNames.size(),0);
        } else {
            totalVariantNumber++;
           // if (totalVariantNumber % reportProgressEvery == 0) reportProgessVCF(totalVariantNumber, startTime);
            
            fields = split(line, '\t');
            VariantInfo v(fields); if (v.onlyIndel) continue; // Only consider SNPs
            
            string offspringGT; bool bInOffspring = false;
            if (offspringPosToGT.count(v.posInt) == 1) {
                offspringGT = offspringPosToGT.at(v.posInt); bInOffspring = true;
                numSharedBetweenVCFs++;
            }
            
            vector<string> genotypes(fields.begin()+NUM_VCF_NON_GENOTYPE_COLUMNS,fields.end());
            // Go through the genotypes and record all SNP alleles
            for (vector<string>::size_type i = 0; i != genotypes.size(); i++) {
                if (genotypes[i][0] == '.' || genotypes[i][2] == '.') continue;
                
                char firstAllele; char secondAllele;
                if (genotypes[i][0] == '0') firstAllele = v.refAllele[0];
                else firstAllele = v.altAlleles[0][0];
                if (genotypes[i][2] == '0') secondAllele = v.refAllele[0];
                else secondAllele = v.altAlleles[0][0];
                
             //   std::cout << "firstAllele: " << firstAllele << std::endl;
             //   std::cout << "secondAllele: " << firstAllele << std::endl;
            
                if (bInOffspring) {
                    char offspringAllele1 = offspringGT[0];
                    char offspringAllele2 = offspringGT[1];
                    if (offspringAllele1 != firstAllele && offspringAllele1 != secondAllele && offspringAllele2 != firstAllele && offspringAllele2 != secondAllele) {
                        numViolations[i]++; continue;
                    }
                }
            }
            genotypes.clear(); genotypes.shrink_to_fit();
        }
    }
    
    print_vector(numViolations, std::cout);
    std::cout << "INFO: DONE. Number of variants shared between the VCFs: " << numSharedBetweenVCFs << std::endl;
    
    return 0;
}



void parseVioOptions(int argc, char** argv) {
    bool die = false; string regionArgString; vector<string> regionArgs;
    vector<string> windowSizeStep;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'n': arg >> opt::runName; break;
            case 'g': opt::useGenotypeProbabilities = true; break;
            case 'h':
                std::cout << VIO_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (argc - optind < 2) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > 2)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << VIO_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::offspringVcfFile = argv[optind++];
    opt::parentsVcfFile = argv[optind++];
}
