#include "analyticalPlacer.hpp"

using namespace std;

void analyticalPlacer::analyticalPlacement(string path) // main algorithm for part 1
{
    processFile(path);
    createNets(blocknets);
    createAdjacencyMatrix();
    matrix_UMF matrixToSolve = createMatrix();
    umfsolver(matrixToSolve);
    HPWL = calculateHPWL();
    cout << HPWL << endl;
}

double analyticalPlacer::calculateHPWL() // calculating hpwl for whole system
{
    double ret = 0;
    for (auto net : nets)
    {

        double lowestX = numeric_limits<double>::max();
        double highestX = numeric_limits<double>::min();

        double lowestY = numeric_limits<double>::max();
        double highestY = numeric_limits<double>::min();
        for (int blockID : net.second)
        {
            if (blockCoordinates[blockID].first < lowestX)
            {
                lowestX = blockCoordinates[blockID].first;
            }
            if (blockCoordinates[blockID].first > highestX)
            {
                highestX = blockCoordinates[blockID].first;
            }
            if (blockCoordinates[blockID].second < lowestY)
            {
                lowestY = blockCoordinates[blockID].second;
            }
            if (blockCoordinates[blockID].second > highestY)
            {
                highestY = blockCoordinates[blockID].second;
            }
        }
        ret += (highestX - lowestX) + (highestY - lowestY);
    }
    return ret;
    // return edge weights
}

map<int, pair<double, double>> analyticalPlacer::umfsolver(matrix_UMF matrixToSolve) // solving the matrix
{
    double *null = (double *)NULL;
    int i;
    int n = nonfixedblockids.size();
    double xid_x[n];
    double xid_y[n];

    void *Symbolic, *Numeric;
    (void)umfpack_di_symbolic(n, n, matrixToSolve.A.Ap.data(), matrixToSolve.A.Ai.data(), matrixToSolve.A.Ax.data(), &Symbolic, null, null);
    (void)umfpack_di_numeric(matrixToSolve.A.Ap.data(), matrixToSolve.A.Ai.data(), matrixToSolve.A.Ax.data(), Symbolic, &Numeric, null, null);

    umfpack_di_free_symbolic(&Symbolic);
    (void)umfpack_di_solve(UMFPACK_A, matrixToSolve.A.Ap.data(), matrixToSolve.A.Ai.data(), matrixToSolve.A.Ax.data(), xid_x, matrixToSolve.Bx.data(), Numeric, null, null);
    (void)umfpack_di_solve(UMFPACK_A, matrixToSolve.A.Ap.data(), matrixToSolve.A.Ai.data(), matrixToSolve.A.Ax.data(), xid_y, matrixToSolve.By.data(), Numeric, null, null);
    umfpack_di_free_numeric(&Numeric);
    for (i = 0; i < n; i++)
    {
        blockCoordinates[nonfixedblockids[i]] = make_pair(xid_x[i], xid_y[i]);
    }
    for (auto i : blockCoordinates)
    {
        cout << i.first << "(" << i.second.first << ", " << i.second.second << ")" << endl;
    }
    return blockCoordinates;
}

void analyticalPlacer::createNets(vector<vector<int>> blocks) // create nets to easily parse between
{
    // i[0] - blockid
    // i[1] - block type
    vector<int> block;
    for (int i = 0; i < blocks.size(); i++)
    {
        block = blocks[i];
        for (int j = 2; j < block.size() - 1; j++)
        {
            if (nets.find(block[j] - 1) != nets.end())
            {
                nets[block[j] - 1].push_back(block[0]);
            }
            else
            {
                nets[block[j] - 1] = {block[0]};
            }
        }
    }
}
void analyticalPlacer::createAdjacencyMatrix() // Create adjacency matrix to keep track of weights between nodes
{
    edgeWeights =
        vector<vector<double>>(blocknets.size() + 1,
                               vector<double>(blocknets.size() + 1, 0));

    int currNetID;
    int pCount;
    for (auto net : nets)
    {
        currNetID = net.first;
        pCount = net.second.size();
        for (int i = 0; i < net.second.size(); i++)
        {
            for (int j = 0; j < net.second.size(); j++)
            {
                if (i != j)
                {
                    edgeWeights[net.second[i]][net.second[j]] += static_cast<double>(2) / pCount;
                }
            }
        }
    }
}
double analyticalPlacer::calculateTotalWeights(int cell, vector<vector<double>> weights)
{
    double tot = 0;
    for (auto val : weights[cell])
    {
        tot += val;
    }
    return tot;
}
matrix_UMF analyticalPlacer::createMatrix()
{
    for (auto net : nets)
    {
        vector<int> ids = net.second;
        for (int id : ids)
        {
            if (fixedblocks.find(id) != fixedblocks.end())
            {
            }
            else
            {
                vector<int>::iterator it = find(nonfixedblockids.begin(), nonfixedblockids.end(), id);
                if (it == nonfixedblockids.end())
                {
                    nonfixedblockids.push_back(id);
                }
            }
        }
    }

    // First create B matrix
    vector<double> Bx = vector<double>(nonfixedblockids.size(), 0);
    vector<double> By = vector<double>(nonfixedblockids.size(), 0);
    vector<vector<double>> A = vector<vector<double>>(nonfixedblockids.size(), vector<double>(nonfixedblockids.size(), 0));
    for (int i = 0; i < nonfixedblockids.size(); i++)
    {
        vector<double> row = edgeWeights[nonfixedblockids[i]];
        for (int colidx = 0; colidx < row.size(); colidx++)
        {
            if (fixedblocks.find(colidx) != fixedblocks.end())
            {
                int connBlock = colidx;
                bool connected = (row[connBlock] > 0) ? true : false;
                if ((row[colidx] > 0))
                {

                    Bx[i] += fixedXCoords[connBlock] * edgeWeights[nonfixedblockids[i]][colidx];

                    By[i] += fixedYCoords[connBlock] * edgeWeights[nonfixedblockids[i]][colidx];
                }
            }
        }
    }

    // Second create A matrix
    // Create actual matrix

    for (int i = 0; i < A.size(); i++)
    {
        for (int j = 0; j < A.size(); j++)
        {
            if (i == j)
            {
                A[i][j] = calculateTotalWeights(nonfixedblockids[i], edgeWeights);
            }
            else
            {
                A[i][j] = (edgeWeights[nonfixedblockids[i]][nonfixedblockids[j]] != 0) ? -edgeWeights[nonfixedblockids[i]][nonfixedblockids[j]] : 0;
            }
        }
    }

    matrix_UMF ret = {matrixToCCForm(A), Bx, By};

    return ret;
}
ACCF analyticalPlacer::matrixToCCForm(vector<vector<double>> matrix) // Convert NxN matrix to compressed column form
{
    vector<int> Ap;
    vector<int> Ai;
    vector<double> Ax;

    int rows = matrix.size();
    int cols = matrix.size();

    Ap.push_back(0);
    int nz = 0;
    for (auto row : matrix)
    {
        for (auto val : row)
        {
            cout << val << " ";
        }
        cout << endl;
    }

    for (int col = 0; col < cols; ++col)
    {
        for (int row = 0; row < rows; ++row)
        {
            if (matrix[row][col] != 0.0)
            {
                Ai.push_back(row);
                Ax.push_back(matrix[row][col]);
                nz++;
            }
        }
        Ap.push_back(nz);
    }
    ACCF ret = {Ap, Ai, Ax};
    return ret;
}

void analyticalPlacer::processFile(const string &path) // Process input file
{
    ifstream infile(path);
    if (!infile)
    {
        cerr << "Error: could not open file!" << endl;
        return;
    }
    string line;
    bool nets = true;
    while (getline(infile, line))
    {
        stringstream ss(line);
        vector<int> values;
        int value;

        while (ss >> value)
        {
            values.push_back(value);
        }

        if (values.size() == 1 && values[0] == -1)
        {
            break;
        }
        blocknets.push_back(values);
        blocktypes[values[0]] = values[1];
    }
    fixedXCoords = vector<int>(blocknets.size() + 1, -1);
    fixedYCoords = vector<int>(blocknets.size() + 1, -1);
    while (getline(infile, line))
    {
        stringstream ss(line);
        vector<int> values;
        int value;

        while (ss >> value)
        {
            values.push_back(value);
        }
        if (values.size() == 1 && values[0] == -1)
        {
            break;
        }
        fixedblocks.insert(values[0]);
        fixedXCoords[values[0]] = values[1];
        fixedYCoords[values[0]] = values[2];
        blockCoordinates[values[0]] = make_pair(values[1] * 1.0, values[2] * 1.0);
    }
}

template <typename T>
void analyticalPlacer::print2Dvector(vector<vector<T>> vec) // For debugging
{
    for (auto i : vec)
    {
        for (auto j : i)
        {
            cout << j << " ";
        }
        cout << endl;
    }
}
//------------------------------------- Part 2 ------------------------------------------

void analyticalPlacer::analyticalPlacementv2(string path, double factor)
{
    processFile(path);
    createNets(blocknets);
    createAdjacencyMatrix();
    adjustFixedWeights(factor);
    matrix_UMF matrixToSolve = createMatrix();
    umfsolver(matrixToSolve);
    HPWL = calculateHPWL();
    cout << HPWL << endl;
}

void analyticalPlacer::adjustFixedWeights(double factor) // increase weight between fixed and movable cells
{
    for (auto block : fixedblocks)
    {
        for (int i = 0; i < edgeWeights.size(); i++)
        {
            edgeWeights[i][block] = edgeWeights[i][block] * factor;
            edgeWeights[block][i] = edgeWeights[block][i] * factor;
        }
    }
}

//------------------------------------- Part 3 ------------------------------------------

double analyticalPlacer::daravspreading() // as described in algorithm 1
{
    bins.clear();
    int ctr = 0;
    splitIntoBins();
    int iter = 2;
    int psi;
    double factor = 1;
    blockCoordinates_new = blockCoordinates;
    blockCoordinates_orig = blockCoordinates;
    overflowed = findOverfilledBins();
    while (!overflowed.empty())
    {
        overflowed = findOverfilledBins();
        // cout << overflowed.size() << endl;
        psi = factor * (iter * iter);
        if (overflowed.size() == 0)
        {
            blockCoordinates = blockCoordinates_new;
            double ret = 0;
            for (int i : nonfixedblockids)
            {
                ret += (fabs(blockCoordinates_orig[i].first - blockCoordinates_new[i].first) + fabs(blockCoordinates_orig[i].second - blockCoordinates_new[i].second));
            }
            return ret;
        }
        sort(overflowed.begin(), overflowed.end(),
             [](const pair<pair<int, int>, int> &a, const pair<pair<int, int>, int> &b)
             { return a.second < b.second; });
        for (pair<pair<int, int>, int> bi : overflowed)
        {

            vector<vector<pair<int, int>>> Pbi = identifyCandidatePaths(bi.first, psi);
            sort(Pbi.begin(), Pbi.end(),
                 [this](const vector<pair<int, int>> &a, const vector<pair<int, int>> &b)
                 { return cost[a] < cost[b]; });
            for (vector<pair<int, int>> pk : Pbi)
            {
                if (supply(bi.first) > 0)
                {
                    cellMove(pk, psi);
                }
            }
        }
        iter++;
    }
    blockCoordinates = blockCoordinates_new;
    double ret = 0;
    for (int i : nonfixedblockids)
    {
        ret += (fabs(blockCoordinates_orig[i].first - blockCoordinates_new[i].first) + fabs(blockCoordinates_orig[i].second - blockCoordinates_new[i].second));
    }
    return ret;
}

vector<vector<pair<int, int>>> analyticalPlacer::identifyCandidatePaths(pair<int, int> bi, int psi) // as described in algorithm 2
{
    int demand = 0;
    vector<vector<pair<int, int>>> ret;
    ret.clear();
    map<pair<int, int>, bool> visited;
    visited.clear();
    for (auto bin : bins)
    {
        visited[bin.first] = false;
    }
    visited[bi] = true;
    vector<pair<int, int>> emptyP;
    emptyP.push_back(bi);
    queue<vector<pair<int, int>>> Q;
    Q.push(emptyP);
    while (!Q.empty() || demand < supply(bi))
    {
        if (Q.empty())
        {
            break;
        }
        vector<pair<int, int>> P = Q.front();
        pair<int, int> tailbin = P[P.size() - 1];
        for (pair<int, int> nDirection : directions)
        {
            pair<int, int> bk = tailbin;
            bk.first += nDirection.first;
            bk.second += nDirection.second;
            if ((0 <= bk.first && bk.first < 30) && (0 <= bk.second && bk.second < 30))
            {
                if (visited[bk])
                {
                    continue;
                }
                int ccost = computecost(tailbin, bk, psi);
                if (ccost < INT_MAX)
                {
                    vector<pair<int, int>> pcopy = P;
                    cost[pcopy] = cost[P] + ccost;
                    pcopy.push_back(bk);
                    visited[bk] = true;
                    if (bins[bk].size() == 0)
                    {
                        ret.push_back(pcopy);
                        demand++;
                    }
                    else
                    {
                        Q.push(pcopy);
                    }
                }
            }
        }
        Q.pop();
    }
    return ret;
}
void analyticalPlacer::cellMove(vector<pair<int, int>> P, int psi) // as described in algorithm 3
{
    int i = 0;
    stack<pair<int, int>> S;
    pair<int, int> Vsink;
    pair<int, int> Vsrc = P[i];
    S.push(P[i]);
    i++;
    while (i != P.size())
    {
        Vsink = P[i];
        int cost = computecost(Vsrc, Vsink, psi);
        if (cost == INT_MAX)
        {
            return;
        }
        Vsrc = Vsink;
        S.push(P[i]);
        i++;
    }
    Vsink = S.top();
    S.pop();
    while (!S.empty())
    {
        Vsrc = S.top();
        S.pop();
        vector<pair<int, double>> blockDistances;
        for (int blockID : bins[Vsrc])
        {
            double quadraticdist = sqrt(pow((Vsink.first - blockCoordinates_orig[blockID].first), 2) + pow((Vsink.second - blockCoordinates_orig[blockID].second), 2));
            blockDistances.push_back({blockID, quadraticdist});
        }
        sort(blockDistances.begin(), blockDistances.end(), [](const pair<int, double> &a, const pair<int, double> &b)
             { return a.second < b.second; });
        if (blockDistances.size() == 0)
        {
            break;
        }
        auto block = blockDistances[0];
        int blockID = block.first;
        pair<double, double> coords = blockCoordinates_new[blockID];
        int xbin = floor(coords.first);
        int ybin = floor(coords.second);
        bins[{xbin, ybin}].erase(blockID);
        if (blocktypes[blockID] == 1)
        {
            if (xbin == 0)
            {
                bins[{xbin + 1, ybin}].erase(blockID);
            }
            else if (xbin < 29)
            {
                int xround = round(coords.first);
                if (xround > xbin)
                {
                    bins[{xbin + 1, ybin}].erase(blockID);
                }
                else
                {
                    bins[{xbin - 1, ybin}].erase(blockID);
                }
            }
            else
            {
                bins[{xbin - 1, ybin}].erase(blockID);
            }
        }

        blockCoordinates_new[block.first] = {blockCoordinates_new[block.first].first - Vsrc.first + Vsink.first, blockCoordinates_new[block.first].second - Vsrc.second + Vsink.second};

        coords = blockCoordinates_new[blockID];
        xbin = floor(coords.first);
        ybin = floor(coords.second);
        bins[{xbin, ybin}].insert(blockID);
        if (blocktypes[blockID] == 1)
        {
            if (xbin == 0)
            {
                bins[{xbin + 1, ybin}].insert(blockID);
            }
            else if (xbin < 29)
            {
                int xround = round(coords.first);
                if (xround > xbin)
                {
                    bins[{xbin + 1, ybin}].insert(blockID);
                }
                else
                {
                    bins[{xbin - 1, ybin}].insert(blockID);
                }
            }
            else
            {
                bins[{xbin - 1, ybin}].insert(blockID);
            }
        }
        blockCoordinates = blockCoordinates_new;
        Vsink = Vsrc;
    }
}
int analyticalPlacer::supply(pair<int, int> bi)
{
    return max(0, usage(bi) - capacity(bi));
}

int analyticalPlacer::usage(pair<int, int> bi)
{ // number of blocks in bin i
    return bins[bi].size();
}
int analyticalPlacer::capacity(pair<int, int> bi)
{
    return 1;
}
double analyticalPlacer::computecost(pair<int, int> tailbin, pair<int, int> bk, int psi) // compute cost function
{
    double ret = INT_MAX;
    for (auto block : bins[tailbin])
    {
        double quadraticdist = sqrt(pow((bk.first - blockCoordinates_orig[block].first), 2) + pow((bk.second - blockCoordinates_orig[block].second), 2));
        if (quadraticdist < psi)
        {
            ret = min(ret, quadraticdist);
        }
    }
    return ret;
}
void analyticalPlacer::splitIntoBins() // Place each block into its respective bin
{
    bins.clear();
    for (int i = 0; i < 30; i++)
    {
        for (int j = 0; j < 30; j++)
        {
            bins[{i, j}] = set<int>();
        }
    }
    for (auto blockID : nonfixedblockids)
    {
        pair<double, double> coords = blockCoordinates_new[blockID];
        int xbin = floor(coords.first);
        int ybin = floor(coords.second);
        bins[{xbin, ybin}].insert(blockID);
        if (blocktypes[blockID] == 1)
        {
            if (xbin == 0)
            {
                bins[{xbin + 1, ybin}].insert(blockID);
            }
            else if (xbin < 29)
            {
                int xround = round(coords.first);
                if (xround > xbin)
                {
                    bins[{xbin + 1, ybin}].insert(blockID);
                }
                else
                {
                    bins[{xbin - 1, ybin}].insert(blockID);
                }
            }
            else
            {
                bins[{xbin - 1, ybin}].insert(blockID);
            }
        }
    }
}
vector<pair<pair<int, int>, int>> analyticalPlacer::findOverfilledBins()
{
    splitIntoBins();
    vector<pair<pair<int, int>, int>> ret;
    ret.clear();
    for (auto bin : bins)
    {
        if (supply(bin.first) > 0)
        {
            ret.push_back({bin.first, bin.second.size()});
        }
    }
    return ret;
}
matrix_UMF analyticalPlacer::createAnchoredMatrix(int weight)
{

    for (auto net : nets)
    {
        vector<int> ids = net.second;
        for (int id : ids)
        {
            if (fixedblocks.find(id) != fixedblocks.end())
            {
            }
            else
            {
                vector<int>::iterator it = find(nonfixedblockids.begin(), nonfixedblockids.end(), id);
                if (it == nonfixedblockids.end())
                {
                    nonfixedblockids.push_back(id);
                }
            }
        }
    }

    // First create B matrix
    vector<double> Bx = vector<double>(nonfixedblockids.size(), 0);
    vector<double> By = vector<double>(nonfixedblockids.size(), 0);
    vector<vector<double>> A = vector<vector<double>>(nonfixedblockids.size(), vector<double>(nonfixedblockids.size(), 0));
    for (int i = 0; i < nonfixedblockids.size(); i++)
    {
        vector<double> row = edgeWeights[nonfixedblockids[i]];
        for (int colidx = 0; colidx < row.size(); colidx++)
        {
            if (fixedblocks.find(colidx) != fixedblocks.end())
            {
                int connBlock = colidx;
                if ((row[colidx] > 0))
                {

                    Bx[i] += fixedXCoords[connBlock] * edgeWeights[nonfixedblockids[i]][colidx];

                    By[i] += fixedYCoords[connBlock] * edgeWeights[nonfixedblockids[i]][colidx];
                }
            }
        }
        Bx[i] += (blockCoordinates_new[nonfixedblockids[i]].first * weight);
        By[i] += (blockCoordinates_new[nonfixedblockids[i]].second * weight);
    }

    // Second create A matrix

    for (int i = 0; i < A.size(); i++)
    {
        for (int j = 0; j < A.size(); j++)
        {
            if (i == j)
            {
                A[i][j] = calculateTotalWeights(nonfixedblockids[i], edgeWeights);
                A[i][j] += (weight);
            }
            else
            {
                A[i][j] = (edgeWeights[nonfixedblockids[i]][nonfixedblockids[j]] != 0) ? -edgeWeights[nonfixedblockids[i]][nonfixedblockids[j]] : 0;
            }
        }
    }

    matrix_UMF ret = {matrixToCCForm(A), Bx, By};

    return ret;
}

void analyticalPlacer::analyticalPlacementAnchored(string path, double factor)
{
    edgeWeights.clear();
    blocktypes.clear();
    nets.clear();
    blockCoordinates.clear();
    blocknets.clear();
    fixedblocks.clear();
    fixedXCoords.clear();
    fixedYCoords.clear();
    processFile(path);
    createNets(blocknets);
    createAdjacencyMatrix();
    adjustFixedWeights(factor);
    matrix_UMF matrixToSolve = createAnchoredMatrix(0);
    umfsolver(matrixToSolve);
    HPWL = calculateHPWL();
    return;
}