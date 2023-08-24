//
//  generalUtils.hpp
//  Hi-reComb
//
//  Created by Milan Malinsky on 02.09.22.
//

#ifndef generalUtils_hpp
#define generalUtils_hpp

#include <stdio.h>
#include <iostream>
#include <map>
#include <vector>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <getopt.h>
#include <assert.h>
#include <algorithm>
#include <random>
#include <ctime>
#include "gzstream.hpp"

using std::string;
#define PROGRAM_BIN "Hi-reComb"
#define PACKAGE_BUGREPORT "millanek@gmail.com"
#define GZIP_EXT ".gz"

#define HelpOption      "       -h, --help                              display this help and exit\n"
#define RunNameOption   "       -n, --run-name                          run-name will be included in the output file name(s)\n"

// VCF format constant
static const int NUM_VCF_NON_GENOTYPE_COLUMNS=9;  // 8 mendatory columns + 1 column with definition of the genotype columns

double stringToDouble(std::string s);
std::vector<std::string> split(const std::string &s, char delim);
std::vector<double> splitToDouble(const std::string &s, char delim);
void splitToDouble(const std::string &s, char delim, std::vector<double> &elems);
void split(const std::string &s, char delim, std::vector<std::string> &elems);
void assertFileOpen(std::ifstream& fh, const std::string& fn);
void assertFileOpen(std::ofstream& fh, const std::string& fn);
std::istream* createReader(const std::string& filename, std::ios_base::openmode mode = std::ios_base::in);
std::ostream* createWriter(const std::string& filename, std::ios_base::openmode mode = std::ios_base::out);

template <class T> double vector_sum(T& vector) {
    double sum = 0;
    for (int i = 0; i < vector.size(); i++) {
        sum += vector[i];
    }
    return sum;
}

template <class T> double vector_average(T& vector) {
    double sum = vector_sum(vector);
    double average = (double)sum / (double)vector.size();
    return average;
}

// Converting numbers (int, double, size_t, and char) to string
template <typename T> std::string numToString(T i) {
    std::string ret;
    std::stringstream out;
    out << i;
    ret = out.str();
    return ret;
}

// Print an arbitrary vector to a file
template <class T> void print_vector(T vector, std::ostream& outFile, char delim = '\t', bool endLine = true) {
    for (int i = 0; i < vector.size(); i++) {
        if (i == (vector.size()-1)) {
            if (endLine) outFile << vector[i] << std::endl;
            else outFile << vector[i];
        } else {
            outFile << vector[i] << delim;
        }
    }
}

class VariantInfo {
public:
    VariantInfo(const std::vector<string>& VCFfields) {
        chr = VCFfields[0]; posInt = atoi(VCFfields[1].c_str());
        refAllele = VCFfields[3]; altAlleles = split(VCFfields[4], ',');
        
        if (refAllele.length() > 1) onlyIndel = true;
        else SNPAlleleIndices.push_back(0);
        
        std::vector<std::string>::iterator it = std::find(altAlleles.begin(), altAlleles.end(), "*");
        if (it != altAlleles.end()) starPos = (int) std::distance(altAlleles.begin(), it);
        else starPos = -1; // There is no star among the alternative alleles
        
        for (int i = 0; i < altAlleles.size(); i++) {
            if (altAlleles[i].length() == 1 && i != starPos) SNPAlleleIndices.push_back(i+1);
        }
        if (SNPAlleleIndices.size() == 0) onlyIndel = true;
        
    }
    
    string chr;
    int posInt;
    string refAllele;
    std::vector<string> altAlleles;
    std::vector<int> SNPAlleleIndices;
    
    bool onlyIndel = false;
    
private:
    int starPos;
};

inline void reportProgessVCF(const int variantsProcessed, const std::clock_t startTime) {
    double durationOverall = ( std::clock() - startTime ) / (double) CLOCKS_PER_SEC;
    std::cout << "Processed " << variantsProcessed << " variants in " << durationOverall << "secs" << std::endl;
}

inline char getAllele(const string& gt, const int alleleNum, const VariantInfo& v) {
    char allele;
    if (gt[alleleNum] == '0') allele = v.refAllele[0];
    else allele = v.altAlleles[0][0];
    return allele;
}

#endif /* generalUtils_hpp */
