#include "lib.h"
using namespace std;
void FMM::DefineGrid(){

    nx = 200;
    ny = 200;
    Lx = 1e2;
    Ly = 1e2;

    dx=Lx/nx;
    dy=Ly/ny;

    numOfNode=(nx+1)*(ny+1);
    numOfElm=nx*ny;

    ifstream ifs_node("node.txt");
    string str1;
    while(getline(ifs_node, str1)){
        istringstream ss(str1);
        string tmp;
        vector<double> tmp_x;
        for(int i=0; i<2; i++){
            getline(ss, tmp, ',');
            tmp_x.push_back(stod(tmp));
        }
        x_vessel.push_back(tmp_x);
    }
    ifs_node.close();

    //for(int i=0; i<x_vessel.size(); i++){
    //    cout << "x = " << x_vessel.at(i).at(0) << "  y = " << x_vessel.at(i).at(1) << endl;
    //}

    Sort_X_Vessel();

    //ofstream outputfile("sort_x.dat");
    //for(int i=0; i<x_vessel.size(); i++){
    //    outputfile << "x = " << x_vessel.at(i).at(0) << "  y = " << x_vessel.at(i).at(1) << endl;
    //}
    //outputfile.close();
    //exit(1);
    double v = 1e0;

    f = 1e0/v;

    string str;
    ifstream ifs("goal.dat");
    int tmp6 = 0;
   // vector<vector<double>> goal;
    while(getline(ifs,str)){
        istringstream ss(str); 
        string tmp5;
        vector<double> tmp_goal;
        goal.emplace_back();
        for(int j=0; j<2; j++){
            getline(ss, tmp5, ' ');
            goal.at(tmp6).push_back(stod(tmp5));
        }
        tmp6++;
    }
    //cout << "goal1のi" << goal.at(0).at(0) << "  goal1のj" << goal.at(0).at(1) << endl;
    //cout << "goal2のi" << goal.at(1).at(0) << "  goal2のj" << goal.at(1).at(1) << endl;
    //exit(1);

    
    x.resize(numOfNode, vector<double>(2, 0));
    T.resize(ny+1, vector<double>(nx+1, 0));
    element.resize(numOfElm);
    T_1D.resize(numOfNode);
    lambda.resize(ny+1, vector<int>(nx+1, 0));

    for(int ic=0;ic<numOfElm;ic++){
        element[ic].node.resize(4);
    }
    
    int tmp3 = 0;
    for(int i=0; i<ny+1; i++){
        for(int j=0; j<nx+1; j++){
            x.at(tmp3).at(0)=j*dx;
            x.at(tmp3).at(1)=i*dy;
            tmp3++;
        }
    }

    tmp3=0;
    for(int i=0;i<ny;i++){
        for(int j=0;j<nx;j++){
        element[tmp3].node[0]=j   +i*(nx+1);
        element[tmp3].node[1]=j+1 +i*(nx+1);
        element[tmp3].node[2]=j+1 +(i+1)*(nx+1);
        element[tmp3].node[3]=j   +(i+1)*(nx+1);
        tmp3++;
        }
    }
    
}

void FMM::Sort_X_Vessel(){
    for(int i=0; i<x_vessel.size(); i++){
        if(x_vessel.at(i).at(1) == 0){
            if(x_vessel.at(0).at(1) == 0){
                if(x_vessel.at(i).at(0) < x_vessel.at(0).at(0)){
                    swap(x_vessel.at(i), x_vessel.at(0));
                }else{
                    continue;
                }
            }else{
                swap(x_vessel.at(i), x_vessel.at(0));
            }
        }
    }

    //double min_distance;
    //double pow_xy = pow(x_vessel.at(0).at(0)-x_vessel.at(1).at(0), 2) + pow(x_vessel.at(0).at(1)-x_vessel.at(1).at(1), 2);
    //min_distance = sqrt(pow_xy);

    for(int i=0; i<x_vessel.size(); i++){
 
      if(x_vessel.at(i).at(1) == 100){

        for(int k=1; k<x_vessel.size(); k++){
            if(x_vessel.at(k).at(1) == 0){
              swap(x_vessel.at(i+1), x_vessel.at(k));
            }
            continue;
        }
        continue;    
      }
        
      double min_distance;
      double pow_xy = pow(x_vessel.at(i).at(0)-x_vessel.at(i+1).at(0), 2) + pow(x_vessel.at(i).at(1)-x_vessel.at(i+1).at(1), 2);
      min_distance = sqrt(pow_xy);


     if(i==x_vessel.size()-2){
        break;
      }

      for(int j=2; j<x_vessel.size()-i; j++){
        double pow_xy_tmp = pow(x_vessel.at(i).at(0)-x_vessel.at(i+j).at(0), 2) + pow(x_vessel.at(i).at(1)-x_vessel.at(i+j).at(1), 2);
          
        double min_distance_tmp = sqrt(pow_xy_tmp);
        if(min_distance_tmp < min_distance){
          min_distance = min_distance_tmp;
          swap(x_vessel.at(i+1), x_vessel.at(i+j));
        }
      }
    }
}

void FMM::CrossedVoxel(){
    for(int i=0; i<element.size(); i++){
          // cout << "1_0" << endl;
        int count = 0;
        bool test;
        test = LeftLineTest(i);
        if(test) count++;
        test = RightLineTest(i);
        if(test) count++;
        test = UpperLineTest(i);
        if(test) count++;
        test = BottomLineTest(i);
        if(test) count++;
        if(count > 0){
            //cout << "「voxel通過」i = " << i << endl;
            InitialSDFNode(i);
        }else{
            //cout << "「voxel通過していない」i = " << i << endl;
            // cout << "1_2" << endl;
            for(int j=0; j<4; j++){
              int mm = element[i].node[j] / (nx+1);
              int nn = element[i].node[j] - mm*(nx+1);
              T.at(mm).at(nn) = 10e2;
              lambda.at(mm).at(nn) = far;
            }
        }
    }

    ofstream outputfile("T_initial.dat");
    for(int i=0; i<nx+1; i++){
        for(int j=0; j<ny+1; j++){
            outputfile << "i= " << i << " j= " << j << " T=" << T.at(i).at(j) << endl;
        }
    }
    outputfile.close();

    int a = 0;
    for(int i=0; i<ny+1; i++){
        for(int j=0; j<nx+1;j++){
            T_1D.at(a) = T.at(j).at(i);
            a++;
        }
    }

    export_vtu("result_initial.vtu");



    for(int i=0; i<ny+1; i++){
        for(int j=0; j<nx+1; j++){
            if(lambda.at(i).at(j) == fix){
                FixGrid(i, j, lambda, T, H);
            }
        }
    }
}

bool FMM::LeftLineTest(int _i){
  for(int i=0; i<x_vessel.size(); i++){
    if(x_vessel.at(i).at(1)==100) continue;
    double x1 = x.at(element[_i].node[0]).at(0);
    double x2 = x.at(element[_i].node[3]).at(0);
    double y1 = x.at(element[_i].node[0]).at(1);
    double y2 = x.at(element[_i].node[3]).at(1);

    double x3 = x_vessel.at(i).at(0);
    double x4 = x_vessel.at(i+1).at(0);
    double y3 = x_vessel.at(i).at(1);
    double y4 = x_vessel.at(i+1).at(1);

    double tc = (x1-x2)*(y3-y1)+(y1-y2)*(x1-x3);
    double td = (x1-x2)*(y4-y1)+(y1-y2)*(x1-x4);

    double ta = (x3-x4)*(y1-y3)+(y3-y4)*(x3-x1);
    double tb = (x3-x4)*(y2-y3)+(y3-y4)*(x3-x2);

    if(tc*td <= 0 && ta*tb <= 0) {
        return true;
      //  cout << "left test passed" << endl;
    }
  }
  return false;
}

bool FMM::RightLineTest(int _i){
  for(int i=0; i<x_vessel.size(); i++){
    if(x_vessel.at(i).at(1)==100) continue;
    double x1 = x.at(element[_i].node[1]).at(0);
    double x2 = x.at(element[_i].node[2]).at(0);
    double y1 = x.at(element[_i].node[1]).at(1);
    double y2 = x.at(element[_i].node[2]).at(1);

    double x3 = x_vessel.at(i).at(0);
    double x4 = x_vessel.at(i+1).at(0);
    double y3 = x_vessel.at(i).at(1);
    double y4 = x_vessel.at(i+1).at(1);

    double tc = (x1-x2)*(y3-y1)+(y1-y2)*(x1-x3);
    double td = (x1-x2)*(y4-y1)+(y1-y2)*(x1-x4);

    double ta = (x3-x4)*(y1-y3)+(y3-y4)*(x3-x1);
    double tb = (x3-x4)*(y2-y3)+(y3-y4)*(x3-x2);

    if(tc*td <= 0 && ta*tb <= 0){
        //cout << "right test passed" << endl;
        return true;  
    } 
    
  }
  return false;
}

bool FMM::UpperLineTest(int _i){
    for(int i=0; i<x_vessel.size(); i++){
      if(x_vessel.at(i).at(1)==100) continue;
      //cout << "_i = " << _i << endl;
      double x1 = x.at(element[_i].node[3]).at(0);
      double x2 = x.at(element[_i].node[2]).at(0);
      double y1 = x.at(element[_i].node[3]).at(1);
      double y2 = x.at(element[_i].node[2]).at(1);

      double x3 = x_vessel.at(i).at(0);
      double x4 = x_vessel.at(i+1).at(0);
      double y3 = x_vessel.at(i).at(1);
      double y4 = x_vessel.at(i+1).at(1);

      double tc = (x1-x2)*(y3-y1)+(y1-y2)*(x1-x3);
      double td = (x1-x2)*(y4-y1)+(y1-y2)*(x1-x4);

      double ta = (x3-x4)*(y1-y3)+(y3-y4)*(x3-x1);
      double tb = (x3-x4)*(y2-y3)+(y3-y4)*(x3-x2);
      if(tc*td <= 0 && ta*tb <= 0){
        //cout << "upper test passed" << endl;
        return true;  
      } 
    }
    return false;
}

bool FMM::BottomLineTest(int _i){
    for(int i=0; i<x_vessel.size(); i++){
      if(x_vessel.at(i).at(1)==100) continue;
      double x1 = x.at(element[_i].node[0]).at(0);
      double x2 = x.at(element[_i].node[1]).at(0);
      double y1 = x.at(element[_i].node[0]).at(1);
      double y2 = x.at(element[_i].node[1]).at(1);

      double x3 = x_vessel.at(i).at(0);
      double x4 = x_vessel.at(i+1).at(0);
      double y3 = x_vessel.at(i).at(1);
      double y4 = x_vessel.at(i+1).at(1);

      double tc = (x1-x2)*(y3-y1)+(y1-y2)*(x1-x3);
      double td = (x1-x2)*(y4-y1)+(y1-y2)*(x1-x4);

      double ta = (x3-x4)*(y1-y3)+(y3-y4)*(x3-x1);
      double tb = (x3-x4)*(y2-y3)+(y3-y4)*(x3-x2);

      if(tc*td <= 0 && ta*tb <= 0){
       // cout << "bottom test passed" << endl;
        return true;  
      } 
    }
    return false;
}


void FMM::InitialSDFNode(int _i){
    for(int j=0; j<4; j++){
         //cout << "1_1_0" << endl;
        double length;
        length = minimumDistance(_i, j);
        int m = element[_i].node[j] / (nx+1);
        int n = element[_i].node[j] - m*(nx+1);
        T.at(m).at(n) = length;
        lambda.at(m).at(n) = fix;
    }


}


double FMM::minimumDistance(int _i, int _j){
    int tmp_point;
    bool is_first = true;
    double min_distance_point;
    double min_distance_point_tmp;

    for(int k=0; k<x_vessel.size(); k++){
         //cout << "1_1_0_0" << endl;

        double pow_tmp = pow(x.at(element[_i].node[_j]).at(0)-x_vessel.at(k).at(0), 2) + pow(x.at(element[_i].node[_j]).at(1)-x_vessel.at(k).at(1), 2);
        min_distance_point_tmp = sqrt(pow_tmp);
        
        if(is_first){
            min_distance_point = min_distance_point_tmp;
            is_first = false;
        }
         if(min_distance_point_tmp < min_distance_point){
            min_distance_point = min_distance_point_tmp;
            tmp_point = k; 
        }
    }
    if(x_vessel.at(tmp_point).at(1) == 0){
        // cout << "1_1_0_1" << endl;
        double bc_x = x_vessel.at(tmp_point+1).at(0) - x_vessel.at(tmp_point).at(0);
        double bc_y = x_vessel.at(tmp_point+1).at(1) - x_vessel.at(tmp_point).at(1);

        double ba_x = x.at(element[_i].node[_j]).at(0) - x_vessel.at(tmp_point).at(0);
        double ba_y = x.at(element[_i].node[_j]).at(1) - x_vessel.at(tmp_point).at(1);

        double pow_bc = pow(bc_x, 2) + pow(bc_y, 2);
        double pow_ba = pow(ba_x, 2) + pow(ba_y, 2);

        double s = (bc_x*ba_x+bc_y*ba_y)/pow_bc;
        double ah = sqrt(pow_ba-s*s*pow_bc);
        if(0 <= s <= 1){
            //if(_i == 2 && _j == 1){
            //    cout << "x" << x.at(element[_i].node[_j]).at(1) << endl;
            //    cout << "x_vessel" <<  x_vessel.at(tmp_point).at(1) << endl;
            //    cout << "ba_x = " << ba_x << endl;
            //    cout << "ba_y = " << ba_y  << endl;
            //    cout << "AB = " << sqrt(pow_ba) << endl;
            //    cout << "BH = " << s*sqrt(pow_bc) << endl;
            //    cout << "s = " << s << endl;
            //    cout << "tmp_point = " << tmp_point << endl;
            //    cout << "_i = " << _i << "_j = " << _j << endl;
            //    cout << "length is nan" << endl;
            //    cout << "ah = " << ah << endl;
            //    exit(1);
            //}
  
            if(isnan(ah)){
                cout << "length is nan" << endl;
                exit(1);
            }
            return ah;  
        } 
    }else if(x_vessel.at(tmp_point).at(1) == 100){
        //cout << "1_1_0_2" << endl;
        double bc_x = x_vessel.at(tmp_point).at(0) - x_vessel.at(tmp_point-1).at(0);
        double bc_y = x_vessel.at(tmp_point).at(1) - x_vessel.at(tmp_point-1).at(1);

        double ba_x = x.at(element[_i].node[_j]).at(0) - x_vessel.at(tmp_point-1).at(0);
        double ba_y = x.at(element[_i].node[_j]).at(1) - x_vessel.at(tmp_point-1).at(1);

        double pow_bc = pow(bc_x, 2) + pow(bc_y, 2);
        double pow_ba = pow(ba_x, 2) + pow(ba_y, 2);

        double s = (bc_x*ba_x+bc_y*ba_y)/pow_bc;
        double ah = sqrt(pow_ba-s*s*pow_bc);
        if(0 <= s <= 1){
            if(isnan(ah)){
                cout << "length is nan" << endl;
                exit(1);
            }
            return ah;  
        } 
    }

    for(int i=0; i<2; i++){
       // cout << "1_1_0_3" << endl;
        double bc_x = x_vessel.at(tmp_point+i).at(0) - x_vessel.at(tmp_point+i-1).at(0);
        double bc_y = x_vessel.at(tmp_point+i).at(1) - x_vessel.at(tmp_point+i-1).at(1);

        double ba_x = x.at(element[_i].node[_j]).at(0) - x_vessel.at(tmp_point+i-1).at(0);
        double ba_y = x.at(element[_i].node[_j]).at(1) - x_vessel.at(tmp_point+i-1).at(1);

        double pow_bc = pow(bc_x, 2) + pow(bc_y, 2);
        double pow_ba = pow(ba_x, 2) + pow(ba_y, 2);

        double s = (bc_x*ba_x+bc_y*ba_y)/pow_bc;
        double ah = sqrt(pow_ba-s*s*pow_bc);
        if(0 <= s <= 1){
            if(isnan(ah)){
                cout << "length is nan" << endl;
                exit(1);
            }
            return ah;  
        } 
    }
    return min_distance_point;
}

bool FMM::IsGoal(vector<vector<int>> &goal,int _i, int _j){
    
    //cout << "is_goal_in" << endl;
    //cout << "goal1のi" << goal.at(0).at(0) << "  goal1のj" << goal.at(0).at(1) << endl;
    //cout << "goal2のi" << goal.at(1).at(0) << "  goal2のj" << goal.at(1).at(1) << endl;
    //cout << _i << _j << endl;
        
    auto is = find_if(
        begin(goal), end(goal),
        [&_i, &_j](const auto& row) {
            //cout << row.at(0) << endl;
            //exit(1);
            return row.at(0)  == _i && row.at(1) == _j;
        }
    );
    if(is == end(goal)){
      return false;
    }else{
      return true;
    }
    /*
    cout << "hello" << endl;
    cout << goal.size() << endl;

    for(int i=0; i<goal.size(); i++){
        if((_i == goal.at(i).at(0)) && (_j == goal.at(i).at(1)) ){
            return true;
        }
    }
    return false;
    */
}





void FMM::FastMarchingMethod(){
   cout << "1" << endl;
    int loop = 0;
    //InitGrid();
    CrossedVoxel();
    cout << "2" << endl;
    bool is_empty_H = false;
    while(!is_empty_H){
        cout << "loop = " << loop << endl;
        cout << "size = " << size << endl;
        DeleteHeap(H, size);
        cout << "3" << endl;
        int i = H.at(size).at(1);
        int j = H.at(size).at(2);
        H.erase(H.begin()+size);
        cout << "i = " << i << " j = " << j << endl;
        FixGrid(i, j, lambda, T, H);
        cout << "size = " << size << endl;
        if(size==0){
            ofstream output_T("T.dat");
            for(int i=0; i<ny+1; i++){
                for(int j=0; j<nx+1; j++){
                    output_T << i << " " << j << " " << T.at(i).at(j) << endl;
                }
            }
            output_T.close();
            is_empty_H = true;
        }
        loop++;
        //if(loop == 3) exit(1);
    }
    //cout << "3" << endl;
}

void FMM::UpdateGrid(int _i, int _j){
    //cout << "update _i = " << _i << endl;
    //cout << "update _j = " << _j  << endl;
    //cout << "2_1_1_1" << endl;
    if(_i == nx && _j != ny && _j != 0){
        //cout << "in1" << endl;
        T_H = T.at(_i-1).at(_j);
        T_V = min(T.at(_i).at(_j-1), T.at(_i).at(_j+1));
    }else if(_j == ny && _i != nx && _i != 0){
        //cout << "in2" << endl;
        T_H = min(T.at(_i-1).at(_j), T.at(_i+1).at(_j));
        T_V = T.at(_i).at(_j-1);
    }else if(_i == 0 && _j != 0 && _j != nx){
        //cout << "in3" << endl;
        T_H = T.at(_i+1).at(_j);
        T_V = min(T.at(_i).at(_j-1), T.at(_i).at(_j+1));   
    }else if(_j == 0 && _i != 0 && _i != nx){
        //cout << "in4" << endl;
        T_H = min(T.at(_i-1).at(_j), T.at(_i+1).at(_j));
        T_V = T.at(_i).at(_j+1);
    }else if(_i == nx && _j == ny){  
       //cout << "in5" << endl;
        T_H = T.at(_i-1).at(_j);
        T_V = T.at(_i).at(_j-1);
    }else if(_i == 0 && _j == 0){
        //cout << "in6" << endl;
        T_H = T.at(_i+1).at(_j);
        T_V = T.at(_i).at(_j+1);
    }else if(_i == nx && _j == 0){
        //cout << "in7" << endl;
        T_H = T.at(_i-1).at(_j);
        T_V = T.at(_i).at(_j+1);
    }else if(_i == 0 && _j == ny){
        T_H = T.at(_i+1).at(_j);
        T_V = T.at(_i).at(_j-1);
    }else{
        T_H = min(T.at(_i-1).at(_j), T.at(_i+1).at(_j));
        T_V = min(T.at(_i).at(_j-1), T.at(_i).at(_j+1));   
    }


    //cout << "2_1_1_2" << endl;
    if( f > abs(T_H-T_V)){
        T.at(_i).at(_j) = (T_H+T_V+sqrt(2*pow(f, 2)-pow((T_H-T_V), 2)))/2;
    }else{
        T.at(_i).at(_j) = f + min(T_H, T_V);
    }
    //cout << "2_1_1_3" << endl;
}

void FMM::FixGrid(int _i, int _j, vector<vector<int>>& _lambda, vector<vector<double>>& _T, vector<vector<int>>& _H){
    _lambda.at(_i).at(_j) = fix;

    //cout << "2_1_00" << endl;
    //cout << "_i = " << _i << "_j = " << _j << endl;
    int step = 0;

    for(int i=_i-1; i<=_i+1; i++){
        for(int j=_j-1; j<=_j+1; j++){
              //cout << "i = " << i << " j = " << j << " step = " << step << endl;
            if(i == nx+1 || j == ny+1 || step == 0 || step == 2 || step == 4 || step == 6 || step == 8 || i < 0 || j < 0){
                //cout << "skip" << endl;
                step++;
                continue;
            }
            //cout << "i = " << i << " j = " << j << " step = " << step << endl;
            bool is_goal = IsGoal(goal, i, j);

            if( !is_goal && _lambda.at(i).at(j) != fix){
                //cout << size << endl;
                //cout << "i = " << i << " j = " << j << " lambda = " << _lambda.at(i).at(j) << endl;
                //cout << "2_1_0" << endl;
                UpdateGrid(i,j);
                //cout << "2_1_1" << endl;
                if(_lambda.at(i).at(j) == far){
                    _lambda.at(i).at(j) = near;
                   // cout << "i = " << i << " j = " << j << " lambda = " << _lambda.at(i).at(j) << endl;
                    //cout << "2_1_2" << endl;
                    
                    H.emplace_back();
                    H.at(size).push_back(tmp);
                    //cout << "2_1_3" << endl;
                    H.at(size).push_back(i);
                    //cout << "2_1_4" << endl;
                    H.at(size).push_back(j);
                    //cout << "2_1_5" << endl;
                    //cout << "size = " << size << "  tmp = " << tmp <<  "  H.at(size).at(0) = " << H.at(size).at(0) << "  H.at(size).at(1) = " <<  H.at(size).at(1)  << "  H.at(size).at(2) = " <<  H.at(size).at(2)  <<  endl;
                    tmp++;
                    size++;
                }else{
                    //cout << "2_1_6" << endl;
                    UpHeap(H,i,j);
                    //cout << "2_1_7" << endl;
                }
            }
            
            step++;
        }
    }
    //cout << "1_1_5" << endl;
}

void FMM::InitGrid(){

    for(int i=0; i<ny+1; i++){
        for(int j=0; j<nx+1; j++){
            bool is_goal2 = IsGoal(goal, i, j);
            if(is_goal2){
                T.at(i).at(j) = 0;
                lambda.at(i).at(j) = fix;
            }else{
                T.at(i).at(j) = 1e10;
                lambda.at(i).at(j) = far;
            }
        }
    }
    //cout << "1_1" << endl;

    for(int i=0; i<ny+1; i++){
        for(int j=0; j<nx+1; j++){
            if(lambda.at(i).at(j) == fix){
                FixGrid(i, j, lambda, T, H);
            }
        }
    }
    //cout << "1_2" << endl;
}
