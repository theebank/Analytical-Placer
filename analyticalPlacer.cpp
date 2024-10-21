#include "analyticalPlacer.hpp"

using namespace std;

int main()
{
    analyticalPlacer ap;
    ap.analyticalPlacement("./tests/cct1.txt");
    cout << "---------------------------------------------------------------------" << endl;
    ap.daravspreading();
    int n = ap.nonfixedblockids.size();
    double xid_x[n];
    double xid_y[n];
    for (auto i : ap.blockCoordinates)
    {
        cout << i.first << "(" << i.second.first << ", " << i.second.second << ")" << endl;
    }
    return 0;
    // print2Dvector(fixedblocks);
}

void analyticalPlacer::analyticalPlacement(string path)
{
    processFile(path);
    createNets(blocknets);
    createAdjacencyMatrix();
    matrix_UMF matrixToSolve = createMatrix();
    umfsolver(matrixToSolve);
    HPWL = calculateHPWL();
    cout << HPWL << endl;
}

double analyticalPlacer::calculateHPWL()
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

map<int, pair<double, double>> analyticalPlacer::umfsolver(matrix_UMF matrixToSolve)
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
    // cout << nonfixedblockids[i] << "(" << xid_x[i] << "," << xid_y[i] << ")" << endl;
}

void analyticalPlacer::createNets(vector<vector<int>> blocks)
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
void analyticalPlacer::createAdjacencyMatrix()
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
    /*
        B Matrix is right portion of equal sign of slide 10 on lecture 5
            [Xa * W(1,A), Xb * W(2,B),...]
        A Matrix is left-most matrix of slide 10 on Lecture 5
            - Matrix is m - by -n, with nz entries
            - represented in this format:
                - int Ap[n+1];
                - int Ai[nz];
                - double Ax[nz];
            - All nonzeros are entries, but an entry may be numerically zero. The row indices of entries in column j are stored in Ai[Ap[j] ... Ap[j+1]-1]. corresponding numerical values are stored in Ax[Ap[j] ... Ax[j+1]-1].
            - No duplicate rππow indices may be present and the row indices in any given column must be sorted in ascending order.
            - First entry Ap[0] must be zero
            - Total number of entries in matrix is thus nz=Ap[n]
            - Except for the fact that extra zero entries may be included there is a unique compressed column representation of any given matrix A
    */
    // Find all blocks that are not fixed
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
                // cout << nonfixedblockids[i] << " " << colidx << " " << edgeWeights[net][nonfixedblockids[i]][colidx] << " " << fixedXCoords[connBlock] << " " << row[colidx] << " " << edgeWeights[net][nonfixedblockids[i]][colidx] * fixedXCoords[connBlock] << endl;
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
ACCF analyticalPlacer::matrixToCCForm(vector<vector<double>> matrix)
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

void analyticalPlacer::processFile(const string &path)
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
void analyticalPlacer::print2Dvector(vector<vector<T>> vec)
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

void analyticalPlacer::adjustFixedWeights(double factor)
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
void analyticalPlacer::daravspreading()
{
    /*
    iter = 0
    while(!B.empty()){
        psi = max allowed movement(iter) //psi is proportional to iter^2
        B = list of overfilled bins in ascending order of supply
            B = identify overflowedBins
            sort bins in order of ascending supply values

        for each bi: B do{
            P(bi) = list of candidate paths for bi
            for each pk: P(bi) do{
                if(supply(bi)>0){
                    move cells along pk
                }
            }
        }
        iter++
    }
    */
    bins.clear();
    overflowed.clear();
    cost.clear();
    int iter = 1;

    int psi;

    overflowed = findOverfilledBins();
    while (!overflowed.empty())
    {
        psi = iter * iter;
        overflowed = findOverfilledBins();
        if (overflowed.size() == 0)
        {
            return;
        }
        for (auto bin : bins)
        {
            if (bin.second > 0)
            {
                overflowed.push_back({{bin.first.first, bin.first.second}, bin.second});
            }
        }
        sort(overflowed.begin(), overflowed.end(),
             [](const pair<pair<int, int>, int> &a, const pair<pair<int, int>, int> &b)
             { return a.second < b.second; });

        for (pair<pair<int, int>, int> bi : overflowed)
        {
            // P(bi) = list of candidate paths for bi
            vector<vector<pair<int, int>>> Pbi = identifyCandidatePaths(bi.first, psi);
            // for each pk: P(bi) do
            for (vector<pair<int, int>> pk : Pbi)
            {
                if (supply(bi.first) > 0)
                {
                    // move cells along pk
                    queue<pair<int, pair<double, double>>> tospread;
                    for (auto block : blockCoordinates)
                    {
                        // if (0 <= (block.second.first - bi.first) <= 1 && 0 <= (block.second.second - bi.second) <= 1)
                        if ((0 <= block.second.first - bi.first.first && block.second.first - bi.first.first <= 1) && (0 <= block.second.second - bi.first.second && block.second.second - bi.first.second <= 1))
                        {
                            if (fixedblocks.find(block.first) != fixedblocks.end())
                            {
                            }
                            else
                            {
                                tospread.push({block.first, {block.second.first - bi.first.first, block.second.second - bi.first.second}});
                            }
                        }
                    }
                    int i = 0;
                    cout << tospread.size() << endl;
                    while (!tospread.empty())
                    {
                        pair<int, pair<double, double>> cell = tospread.front();
                        cout << cell.second.first + pk[i].first << "," << cell.second.second + pk[i].second << "," << cell.first << endl;
                        blockCoordinates[cell.first] = {cell.second.first + pk[i].first, cell.second.second + pk[i].second};
                        bins[{pk[i].first, pk[i].second}] += 1;
                        bins[bi.first] -= 1;
                        i++;
                        tospread.pop();
                    }
                }
            }
        }
        iter++;
    }
}
vector<vector<pair<int, int>>> analyticalPlacer::identifyCandidatePaths(pair<int, int> bi, int psi)
{
    vector<vector<pair<int, int>>> ret;
    map<pair<int, int>, bool> visited;
    //    demand = 0
    int demand = 0;
    // mark all bins as unvisited
    for (auto bin : bins)
    {
        visited[bin.first] = false;
    }
    // visited (bi) = true
    visited[bi] = true;
    // insert bi into empty path P
    vector<pair<int, int>> emptyP;
    emptyP.push_back(bi);
    vector<pair<int, int>> p;
    // add path p to FIFO q
    queue<vector<pair<int, int>>> fifo;
    fifo.push(emptyP);

    pair<int, int> bk;
    // while(!Q.empty() || demand>=supply(bi)){
    while (!fifo.empty() || demand >= supply(bi))
    {
        if (fifo.empty())
        {
            break;
        }
        //     pop p from Q
        p = fifo.front();

        //     tail bin = current end of path p
        pair<int, int> tail = p[p.size() - 1];
        //     for neighbour bins bk of tailbin{
        bk = tail;
        for (pair<int, int> nDirection : directions)
        {

            bk.first += nDirection.first;
            bk.second += nDirection.second;
            if ((0 <= bk.first && bk.first <= 29) && (0 <= bk.second && bk.second <= 29))
            {
                // if bk is visited{
                if (visited[bk])
                {
                    // continue/ignore
                    continue;
                }
                // cout << bk.first << bk.second << endl;

                // ccost = compute cost(tailbin, bk) <- cost of extending path by one bin
                int ccost = computecost(tail, bk, psi);
                // if(cost < inf){
                if (ccost < INT_MAX)
                {
                    // pcopy = a copy of p
                    vector<pair<int, int>> pcopy = p;
                    // cost(pcopy) = cost(p) + ccost
                    cost[pcopy] = cost[p] + ccost;
                    // insert bk to pcopy //new tail
                    pcopy.push_back(bk);
                    // visited(bk) = true
                    visited[bk] = true;
                    // if(bk.empty()){
                    if (bins[bk] == 0)
                    {

                        // insert pcopy to P(bi)
                        ret.push_back(pcopy);
                        // demand++
                        demand++;
                        //}else{
                    }
                    else
                    {
                        // add pcopy into Q
                        fifo.push(pcopy);
                    }
                }
            }
        }
        fifo.pop();
    }

    return ret;
}
int analyticalPlacer::supply(pair<int, int> bi)
{
    int ret = 0;
    ret = max(0, usage(bi) - capacity(bi));
    return ret;
}

int analyticalPlacer::usage(pair<int, int> bi)
{ // number of blocks in bin i
    return bins[bi];
}
int analyticalPlacer::capacity(pair<int, int> bi)
{
    return 1;
}
int analyticalPlacer::computecost(pair<int, int> tailbin, pair<int, int> bk, int psi)
{
    double quadraticdist = pow((bk.first - tailbin.first), 2) + pow((bk.second - tailbin.second), 2);

    if (static_cast<int>(quadraticdist) > psi)
    {
        return INT_MAX;
    }
    else
    {
        // cout << "COST IS: " << static_cast<int>(quadraticdist) << endl;
        return static_cast<int>(quadraticdist);
    }
}
void analyticalPlacer::splitIntoBins()
{
    for (int i = 0; i < 30; i++)
    {
        for (int j = 0; j < 30; j++)
        {
            bins[{i, j}] = 0;
        }
    }
    for (auto blockid : nonfixedblockids)
    {
        pair<double, double> coords = blockCoordinates[blockid];
        int xbin = static_cast<int>(coords.first);
        int ybin = static_cast<int>(coords.second);

        bins[{xbin, ybin}] += 1;
    }
    cout << endl;
}
vector<pair<pair<int, int>, int>> analyticalPlacer::findOverfilledBins()
{
    splitIntoBins();
    vector<pair<pair<int, int>, int>> ret;
    for (auto bin : bins)
    {
        if (bin.second > 1)
        {
            ret.push_back({bin.first, bin.second});
        }
    }
    return ret;
}
