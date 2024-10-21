#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <utility>
#include <map>
#include <limits>
#include <algorithm>
#include "umfpack.h"

#ifndef ANALYTICALPLACER_HPP
#define ANALYTICALPLACER_HPP

using namespace std;

struct ACCF
{
    vector<int> Ap;
    vector<int> Ai;
    vector<double> Ax;
};

struct matrix_UMF
{
    ACCF A;
    vector<double> Bx;
    vector<double> By;
};

class analyticalPlacer
{
public:
    vector<vector<int>> blocknets;
    set<int> fixedblocks;
    vector<int> fixedXCoords;
    vector<int> fixedYCoords;

    vector<int> nonfixedblockids;

    // Adjacency matrix for edge weights
    vector<vector<double>> edgeWeights;
    vector<vector<int>> blockLocations;

    map<int, int> blocktypes;
    map<int, vector<int>> nets;
    map<int, pair<double, double>> blockCoordinates;
    double HPWL;

    void processFile(const string &path);
    template <typename T>
    void print2Dvector(vector<vector<T>> vec);

    double calculateTotalWeights(int cell, vector<vector<double>> weights);

    void createNets(vector<vector<int>> blocks);
    void createAdjacencyMatrix();
    matrix_UMF createMatrix();
    ACCF matrixToCCForm(vector<vector<double>> matrix);
    double calculateHPWL();
    map<int, pair<double, double>> umfsolver(matrix_UMF matrixToSolve);
    void analyticalPlacement(string path);

    void adjustFixedWeights(double factor);
    void analyticalPlacementv2(string path, double factor);

    void daravspreading();
    void identifyCandidatePaths();
};
#endif // ANALYTICALPLACER_HPP