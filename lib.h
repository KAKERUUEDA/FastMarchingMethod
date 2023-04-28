#include <bits/stdc++.h>
#include <vector>
using namespace std;


enum lambda_type{
    fix = 1,
    near = 2,
    far = 3, 
};

class ElementType{
    public:
    vector<int> node;
};

class FMM{
    public:
    int nx, ny, goal_i, goal_j, tmp, tmp2, size, numOfNode, numOfElm;
    double Lx, Ly;
    double dx, dy;
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

    void Sort_X_Vessel();
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
    //void AddHeap(vector<int> &H, int size);
    void UpHeap(vector<vector<int>>& H, int _i, int _j);
    void DeleteHeap(vector<vector<int>>& H, int &size);

    void export_vtu(const string &file);
};