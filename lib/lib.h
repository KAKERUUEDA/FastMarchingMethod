#include <iostream>
#include <vector>
#include <string>
#include <fstream>  
#include <sstream> 
#include <cmath>
#include <algorithm>
using namespace std;

enum lambda_type
{
    fix = 1,
    near = 2,
    far = 3, 
};

class ElementType
{
    public:
    vector<int> node;
};

class FMM
{
public:
    FMM() : tmp(0), tmp2(0), size(0) {}
    ~FMM() {}

    int nx, ny;
    double Lx, Ly;
    double dx, dy;
    int numOfNode, numOfElm;
    int tmp, tmp2;
    int size;
    int loop;
    int goal_i, goal_j;
    double f;
    double T_H, T_V;
    static const int N;

    vector<vector<int>> H;
    vector<double> T_1D;
    vector<vector<double>> x;
    vector<vector<double>> x_vessel;
    vector<ElementType> element;
    vector<vector<double>> T;
    vector<vector<int>> lambda;
    vector<vector<int>> goal;

    void Sort();
    void CrossedVoxel();
    bool LeftLineTest(int _i);
    bool RightLineTest(int _i);
    bool UpperLineTest(int _i);
    bool BottomLineTest(int _i);
    void InitialSDFNode(int _i);
    double minimumDistance(int _i, int _j);

    bool IsGoal(vector<vector<int>> &goal, int _i, int _j);

    void DefineGrid();
    void InitGrid();
    void FixGrid(int _i, int _j, vector<vector<int>>& _lambda, vector<vector<double>>& _T, vector<vector<int>>& _H);
    void UpdateGrid(int _i, int _j);
    void FastMarchingMethod();
    void UpHeap(vector<vector<int>>& H, int _i, int _j);
    void DeleteHeap(vector<vector<int>>& H, int &size);
    void updateT1D();
    void export_vtu(const string &file);
};

