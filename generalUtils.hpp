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

using std::string;
#define PROGRAM_BIN "Hi-reComb"
#define PACKAGE_BUGREPORT "mm21@sanger.ac.uk"
#define GZIP_EXT ".gz"

double stringToDouble(std::string s);
std::vector<std::string> split(const std::string &s, char delim);
std::vector<double> splitToDouble(const std::string &s, char delim);
void splitToDouble(const std::string &s, char delim, std::vector<double> &elems);
void split(const std::string &s, char delim, std::vector<std::string> &elems);
void assertFileOpen(std::ifstream& fh, const std::string& fn);
void assertFileOpen(std::ofstream& fh, const std::string& fn);

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

#endif /* generalUtils_hpp */
