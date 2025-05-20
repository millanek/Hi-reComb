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
#include <cmath>
#include "gzstream.hpp"

using std::string;
using std::vector;
using std::map;

#define PROGRAM_BIN "Hi-reComb"
#define PACKAGE_BUGREPORT "millanek@gmail.com"
#define GZIP_EXT ".gz"

#define HelpOption      "       -h, --help                              display this help and exit\n"
#define RunNameOption   "       -n, --run-name                          run-name will be included in the output file name(s)\n"

// VCF format constant
static const int NUM_VCF_NON_GENOTYPE_COLUMNS=9;  // 8 mendatory columns + 1 column with definition of the genotype columns

double stringToDouble(string s);
vector<string> split(const string &s, char delim);
vector<double> splitToDouble(const string &s, char delim);
void splitToDouble(const string &s, char delim, vector<double> &elems);
void split(const string &s, char delim, vector<string> &elems);
void assertFileOpen(std::ifstream& fh, const string& fn);
void assertFileOpen(std::ofstream& fh, const string& fn);
std::istream* createReader(const string& filename, std::ios_base::openmode mode = std::ios_base::in);
std::ostream* createWriter(const string& filename, std::ios_base::openmode mode = std::ios_base::out);

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
template <typename T> string numToString(T i) {
    string ret;
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
    VariantInfo(const vector<string>& VCFfields) {
        chr = VCFfields[0]; posInt = atoi(VCFfields[1].c_str());
        refAllele = VCFfields[3]; altAlleles = split(VCFfields[4], ',');
        
        if (refAllele.length() > 1) onlyIndel = true;
        else SNPAlleleIndices.push_back(0);
        
        vector<string>::iterator it = std::find(altAlleles.begin(), altAlleles.end(), "*");
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
    vector<string> altAlleles;
    vector<int> SNPAlleleIndices;
    
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


// -------------------------------------   BASIC RECOMBINATION MAP DEFINITION CLASSES  ----------------------------------------
// These are extended in specific use-cases 

class RecombIntervalBase {
public:
    RecombIntervalBase() {};
    
    RecombIntervalBase(int leftIndex, int leftCoordIn, int rightCoordIn, double perBpRateIn) {
        j = leftIndex;
        recombFractionPerBp = perBpRateIn;
        leftCoord = leftCoordIn; rightCoord = rightCoordIn;
        dj = rightCoord - leftCoord + 1;
        rj = recombFractionPerBp * dj;
    };
    
    int j;
    int leftCoord;
    int rightCoord;
    double recombFractionPerBp;
    int dj; double rj;
};


class RecombMapBase {
public:
    
    RecombMapBase() {
        physicalWindowStartEnd.resize(2);
    };
    
    double meanRate;
    vector<double> cummulativeRates; // Cummulative rates in cM
    double mapLength; // Map length in cM
    long long int mapPhysicalLength; // Map length in bp
    int mapPhysicalStart; // Map start in bp position on chromosome
    int mapPhysicalEnd; // Map end in bp position on chromosome
    double meanEffectiveCoverage;
    vector<double> physicalWindowR;
    vector<vector<int>> physicalWindowStartEnd;
    vector<vector<double>> physicalWindowBootstraps;
    
    double getAverageRateForPhysicalWindow(const int start, const int end) {
        
        if (intervalCoordsVectors.empty()) {
            std::cerr << "ERROR: problem calculating average rate for a particular genomic window starting at: " << start << "bp" << std::endl;
            std::cerr << "This is most likely a bug. Please report. " << std::endl;
            exit(EXIT_FAILURE);
        }
        
        // Binary search to find the first interval whose end coordinate is greater
        // or equal to the start of the region in question
        vector<int>::iterator itStart = lower_bound(intervalCoordsVectors[1].begin(),intervalCoordsVectors[1].end(),start);
        int numBPtotal = 0; int numBPthisInterval = 0;
        double sumPerBPvalue = 0; double meanValue = NAN;
        
        if (itStart != intervalCoordsVectors[1].end()) {  // if (start < f[1])    ---  excludind case 1)
            vector<int>::size_type index = std::distance(intervalCoordsVectors[1].begin(), itStart);
            // Sum the lengths
            while (intervalCoordsVectors[0][index] <= end && index < intervalCoordsVectors[0].size()) { // if (f[0] >= end)   ---    excluding case 2)
                double valueThisInterval = intervalPerBPrVector[index];

              //  std::cerr << "valuesThisFeature[0]\t" << valuesThisFeature[0] << std::endl;
                if (intervalCoordsVectors[0][index] < start && intervalCoordsVectors[1][index] <= end)
                    numBPthisInterval = (intervalCoordsVectors[1][index] - start) + 1;
                else if (intervalCoordsVectors[0][index] >= start && intervalCoordsVectors[1][index] <= end)
                    numBPthisInterval = (intervalCoordsVectors[1][index] - intervalCoordsVectors[0][index]);
                else if (intervalCoordsVectors[0][index] >= start && intervalCoordsVectors[1][index] > end)
                    numBPthisInterval = (end - intervalCoordsVectors[0][index]);
                else if (intervalCoordsVectors[0][index] < start && intervalCoordsVectors[1][index] > end)
                    numBPthisInterval = (end - start) + 1;
                
                numBPtotal += numBPthisInterval;
                sumPerBPvalue += valueThisInterval * numBPthisInterval;
                index++;
            }
            
                meanValue = (double)sumPerBPvalue/numBPtotal;
            //std::cerr << "meanValue: " << meanValue << std::endl;
        }
        return meanValue;
    }
    
protected:
    vector<vector<int>> intervalCoordsVectors; vector<double> intervalPerBPrVector;
};



// -------------------------------------    BASIC MATH/STATS  ----------------------------------------

// factorial(x): (x! for non-negative integer x) is defined to be gamma(x+1) (as in R)
inline double factorial(double num) {
    if (num < 0) {
        std::cerr << "Can't compute factorial of a negative number " << num << std::endl;
        exit(1);
    }
    return tgamma(num+1);
}

// Calculates the binomial coefficient (n choose k)
inline int choose(int n, int k) {
    double dResult = factorial(n)/(factorial(k)*factorial(n-k));
    int iResult = (int)round(dResult);
    return iResult;
}

inline long double binomialPMF(const int k, const int n, const double p) {
    return choose(n,k) * pow(p,k) * pow(1  - p, n - k);
}



#endif /* generalUtils_hpp */
