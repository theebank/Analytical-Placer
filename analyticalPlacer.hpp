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
#include <queue>
#include <cmath>
#include <stack>
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
    vector<vector<pair<int, int>>> identifyCandidatePaths(pair<int, int> bi, int psi);
    void cellMove(vector<pair<int, int>> P, int psi);
    int supply(pair<int, int> bi);
    int usage(pair<int, int> bi);
    int capacity(pair<int, int> bi);
    double computecost(pair<int, int> tailbin, pair<int, int> bk, int psi);
    void splitIntoBins();
    vector<pair<pair<int, int>, int>> findOverfilledBins();
    map<pair<int, int>, set<int>> bins; // key = coordinate, value = blocks inside

    map<int, pair<double, double>> blockCoordinates_orig;
    map<int, pair<double, double>> blockCoordinates_new;
    vector<pair<pair<int, int>, int>> overflowed; // key = bin_num, value = supply
    // map<pair<int, int>, int>
    //     bins;                                     // key = coordinate, value = supply
    map<vector<pair<int, int>>, int> cost;
    vector<pair<int, int>> directions = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
};
#endif // ANALYTICALPLACER_HPP