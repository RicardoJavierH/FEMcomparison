//
// Created by victor on 02/09/2020.
//

#ifndef FEMCOMPARISON_INPUTTREATMENT_H
#define FEMCOMPARISON_INPUTTREATMENT_H

#include "DataStructure.h"

void DataInitialization(int argc, char *argv[],PreConfig &pConfig);

void EvaluateEntry(int argc, char *argv[],PreConfig &eData);
void ReadEntry(ProblemConfig &config, PreConfig &preConfig);
void InitializeOutstream(PreConfig &eData,char *argv[]);
void IsInteger(char *argv);
void IsInteger(std::string str);
void isFloat( std::string myString );

void Configure(ProblemConfig &config,int ndiv,PreConfig &pConfig,char *argv[]);
void CharReplace(std::string &str, char find, char replace );
void InitializeSpeedUp(PreConfig &pConfig);
void DirectoryCreation(std::string &resultsFile, std::string &time);
void InitializeOfsSim(PreConfig &pConfig);

template <typename Out>
void split(const std::string &s, char delim, Out result);
std::vector<std::string> split(const std::string &s, char delim);
void ManageOfstream(PreConfig &pConfig, int nThreads, int iterNum);
void ManageOfstream(PreConfig &pConfig);

void OfstreamPath(PreConfig &pConfig);
void ReadSimFile(std::vector<std::string> &cellVec, PreConfig &pConfig);

#endif //FEMCOMPARISON_INPUTTREATMENT_H
