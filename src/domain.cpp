#include "../lib/lib.h"

void FMM::DefineGrid()
{
    nx = 200;
    ny = 200;
    Lx = 1e2;
    Ly = 1e2;

    dx=Lx/nx;
    dy=Ly/ny;

    numOfNode = (nx + 1) * (ny + 1);
    numOfElm = nx * ny;

    ifstream ifs_node("../../input/node.txt");
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

    Sort();

    double v = 1e0;

    f = 1e0/v;

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

void FMM::Sort()
{
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

void FMM::CrossedVoxel()
{
    for(int i=0; i<element.size(); i++){
        int count = 0;
        bool test;
        test = LeftLineTest(i);
        if(test){
            count++;
        } 
        test = RightLineTest(i);
           if(test){
            count++;
        } 
        test = UpperLineTest(i);
           if(test){
            count++;
        } 
        test = BottomLineTest(i);
           if(test){
            count++;
        } 
        if(count > 0){
            InitialSDFNode(i);
        }else{
            for(int j=0; j<4; j++){
              int mm = element[i].node[j] / (nx+1);
              int nn = element[i].node[j] - mm*(nx+1);
              if(T.at(mm).at(nn) != 0e0) continue;
              T.at(mm).at(nn) = 1e4;
              lambda.at(mm).at(nn) = far;
            }
        }
    }

    int a = 0;
    for(int i=0; i<ny+1; i++){
        for(int j=0; j<nx+1;j++){
            T_1D.at(a) = T.at(i).at(j);
            a++;
        }
    }

    for(int i=0; i<ny+1; i++){
        for(int j=0; j<nx+1; j++){
            if(lambda.at(i).at(j) == fix){
                FixGrid(i, j, lambda, T, H);
            }
        }
    }
}

bool FMM::LeftLineTest(int _i)
{
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
    }
  }
  return false;
}

bool FMM::RightLineTest(int _i)
{
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
        return true;  
    } 
    
  }
  return false;
}

bool FMM::UpperLineTest(int _i)
{
    for(int i=0; i<x_vessel.size(); i++){
      if(x_vessel.at(i).at(1)==100) continue;
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
        return true;  
      } 
    }
    return false;
}


void FMM::InitialSDFNode(int _i)
{
    for(int j=0; j<4; j++){
        double length;
        length = minimumDistance(_i, j);
        int m = element[_i].node[j] / (nx+1);
        int n = element[_i].node[j] - m*(nx+1);
        T.at(m).at(n) = length;
        lambda.at(m).at(n) = fix;
    }
} 



double FMM::minimumDistance(int _i, int _j)
{
    int tmp_point;
    double min_distance_point;
    double min_distance_point_tmp;

    double pow_tmp = pow(x.at(element[_i].node[_j]).at(0)-x_vessel.at(0).at(0), 2) + pow(x.at(element[_i].node[_j]).at(1)-x_vessel.at(0).at(1), 2);
    min_distance_point = sqrt(pow_tmp);
    tmp_point = 0;

    for(int k=0; k<x_vessel.size(); k++){
        double pow_tmp = pow(x.at(element[_i].node[_j]).at(0)-x_vessel.at(k).at(0), 2) + pow(x.at(element[_i].node[_j]).at(1)-x_vessel.at(k).at(1), 2);
        min_distance_point_tmp = sqrt(pow_tmp);

        if(min_distance_point_tmp < min_distance_point){
            min_distance_point = min_distance_point_tmp;
            tmp_point = k; 
        }
    }
    if(x_vessel.at(tmp_point).at(1) == 0){
        double bc_x = x_vessel.at(tmp_point+1).at(0) - x_vessel.at(tmp_point).at(0);
        double bc_y = x_vessel.at(tmp_point+1).at(1) - x_vessel.at(tmp_point).at(1);

        double ba_x = x.at(element[_i].node[_j]).at(0) - x_vessel.at(tmp_point).at(0);
        double ba_y = x.at(element[_i].node[_j]).at(1) - x_vessel.at(tmp_point).at(1);

        double pow_bc = pow(bc_x, 2) + pow(bc_y, 2);
        double pow_ba = pow(ba_x, 2) + pow(ba_y, 2);

        double s = (bc_x*ba_x+bc_y*ba_y)/pow_bc;
        double ah = sqrt(pow_ba-s*s*pow_bc);

        if(s >= 0e0 && 1e0 >= s){
            return ah;  
        } 
    }else if(x_vessel.at(tmp_point).at(1) == 100){
        double bc_x = x_vessel.at(tmp_point).at(0) - x_vessel.at(tmp_point-1).at(0);
        double bc_y = x_vessel.at(tmp_point).at(1) - x_vessel.at(tmp_point-1).at(1);

        double ba_x = x.at(element[_i].node[_j]).at(0) - x_vessel.at(tmp_point-1).at(0);
        double ba_y = x.at(element[_i].node[_j]).at(1) - x_vessel.at(tmp_point-1).at(1);

        double pow_bc = pow(bc_x, 2) + pow(bc_y, 2);
        double pow_ba = pow(ba_x, 2) + pow(ba_y, 2);

        double s = (bc_x*ba_x+bc_y*ba_y)/pow_bc;
        double ah = sqrt(pow_ba-s*s*pow_bc);
        if(s >= 0e0 && 1e0 >= s){
            if(isnan(ah)){
                cout << "length is nan" << endl;
                exit(1);
            }
            return ah;  
        } 
    }else if(x_vessel.at(tmp_point).at(1) != 0  && x_vessel.at(tmp_point).at(1) != 100){
        for(int i=0; i<2; i++){
          double bc_x = x_vessel.at(tmp_point+i).at(0) - x_vessel.at(tmp_point+i-1).at(0);
          double bc_y = x_vessel.at(tmp_point+i).at(1) - x_vessel.at(tmp_point+i-1).at(1);

          double ba_x = x.at(element[_i].node[_j]).at(0) - x_vessel.at(tmp_point+i-1).at(0);
          double ba_y = x.at(element[_i].node[_j]).at(1) - x_vessel.at(tmp_point+i-1).at(1);

          double pow_bc = pow(bc_x, 2) + pow(bc_y, 2);
          double pow_ba = pow(ba_x, 2) + pow(ba_y, 2);

          double s = (bc_x*ba_x+bc_y*ba_y)/pow_bc;
          double ah = sqrt(pow_ba-s*s*pow_bc);
          if(s >= 0e0 && 1e0 >= s){
            if(isnan(ah)){
                cout << "length is nan" << endl;
                exit(1);
            }
          return ah;  
          } 
        }
    }

    return min_distance_point;
}

bool FMM::IsGoal(vector<vector<int>> &goal,int _i, int _j)
{
    auto is = find_if(
        begin(goal), end(goal),
        [&_i, &_j](const auto& row) {
            return row.at(0)  == _i && row.at(1) == _j;
        }
    );
    if(is == end(goal)){
      return false;
    }else{
      return true;
    }
}





void FMM::FastMarchingMethod()
{
    loop = 0;
    CrossedVoxel();
    bool is_empty_H = false;

    while(!is_empty_H){
        cout << "loop = " << loop << endl;
        DeleteHeap(H, size);
        int i = H.at(size).at(1);
        int j = H.at(size).at(2);
        H.erase(H.begin()+size);
        FixGrid(i, j, lambda, T, H);
        if(size == 0){
            is_empty_H = true;
        }
        int a = 0;
        for(int i=0; i<ny+1; i++){
          for(int j=0; j<nx+1;j++){
            T_1D.at(a) = T.at(i).at(j);
            a++;
          }
        }

        loop++;
    }
}

void FMM::UpdateGrid(int _i, int _j){
    if(_i == nx && _j != ny && _j != 0){
        T_H = T.at(_i-1).at(_j);
        T_V = min(T.at(_i).at(_j-1), T.at(_i).at(_j+1));
    }else if(_j == ny && _i != nx && _i != 0){
        T_H = min(T.at(_i-1).at(_j), T.at(_i+1).at(_j));
        T_V = T.at(_i).at(_j-1);
    }else if(_i == 0 && _j != 0 && _j != nx){
        T_H = T.at(_i+1).at(_j);
        T_V = min(T.at(_i).at(_j-1), T.at(_i).at(_j+1));   
    }else if(_j == 0 && _i != 0 && _i != nx){
        T_H = min(T.at(_i-1).at(_j), T.at(_i+1).at(_j));
        T_V = T.at(_i).at(_j+1);
    }else if(_i == nx && _j == ny){  
        T_H = T.at(_i-1).at(_j);
        T_V = T.at(_i).at(_j-1);
    }else if(_i == 0 && _j == 0){
        T_H = T.at(_i+1).at(_j);
        T_V = T.at(_i).at(_j+1);
    }else if(_i == nx && _j == 0){
        T_H = T.at(_i-1).at(_j);
        T_V = T.at(_i).at(_j+1);
    }else if(_i == 0 && _j == ny){
        T_H = T.at(_i+1).at(_j);
        T_V = T.at(_i).at(_j-1);
    }else{
        T_H = min(T.at(_i-1).at(_j), T.at(_i+1).at(_j));
        T_V = min(T.at(_i).at(_j-1), T.at(_i).at(_j+1));  
    }

    if( f > fabs(T_H-T_V)){
        T.at(_i).at(_j) = (T_H+T_V+sqrt(2*pow(f, 2)-pow((T_H-T_V), 2)))/2;
    }else{
        T.at(_i).at(_j) = f + min(T_H, T_V);
    }
}

void FMM::FixGrid(int _i, int _j, vector<vector<int>>& _lambda, vector<vector<double>>& _T, vector<vector<int>>& _H)
{
    _lambda.at(_i).at(_j) = fix;
    UpdateGrid(_i,_j);

    int step = 0;

   for(int i=_i-1; i<=_i+1; i++){
       for(int j=_j-1; j<=_j+1; j++){
            if(i == nx+1 || j == ny+1 || step == 0 || step == 2 || step == 4 || step == 6 || step == 8 || i < 0 || j < 0){
                step++;
                continue;
            }

            if( _lambda.at(i).at(j) != fix){
                UpdateGrid(i,j);

                if(_lambda.at(i).at(j) == far){
                    _lambda.at(i).at(j) = near;
                    
                    H.emplace_back();

                    H.at(size).push_back(tmp);
                    H.at(size).push_back(i);
                    H.at(size).push_back(j);

                    tmp++;
                    size++;
                }else{
                    UpHeap(H,i,j);
                }
            }
            
            step++;
        }
    }
}

void FMM::InitGrid()
{
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

    for(int i=0; i<ny+1; i++){
        for(int j=0; j<nx+1; j++){
            if(lambda.at(i).at(j) == fix){
                FixGrid(i, j, lambda, T, H);
            }
        }
    }
}

void FMM::updateT1D()
{
    int a = 0;
    for(int i=0; i<ny+1; i++){
        for(int j=0; j<nx+1;j++){
            T_1D.at(a) = T.at(i).at(j);
            a++;
        }
    }
}

void FMM::export_vtu(const string &file)
{
  FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", numOfNode, numOfElm);
  fprintf(fp, "<Points>\n");
  int offset = 0;
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ic=0;ic<numOfNode;ic++){
    fprintf(fp,"%e %e 0e0\n",x.at(ic).at(0),x.at(ic).at(1));
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Points>\n");

  fprintf(fp, "<Cells>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
  for (int i = 0; i < numOfElm; i++){
    for (int j = 0; j < element[i].node.size(); j++) fprintf(fp, "%d ", element[i].node[j]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
  int num = 0;
  for (int i = 0; i < numOfElm; i++)
  {
    num += element[i].node.size();
    fprintf(fp, "%d\n", num);
  }
  
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (int i = 0; i < numOfElm; i++) fprintf(fp, "%d\n", 5);
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Cells>\n");
  

  fprintf(fp, "<PointData>\n");
   fprintf(fp, "<DataArray type=\"Float64\" Name=\"T_1D\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    for(int ic=0;ic<numOfNode;ic++){
      fprintf(fp,"%e\n",T_1D.at(ic));
    }
    fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</PointData>\n");

  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}
