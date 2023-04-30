#include <iostream>
#include "lib.h"
using namespace std;

void FMM::export_vtu(const std::string &file)
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


int main(){
    FMM fmm;
    fmm.tmp=0;
    fmm.tmp2=0;
    fmm.size=0;
    fmm.DefineGrid();
    fmm.FastMarchingMethod();
    
    int a = 0;
    for(int i=0; i<fmm.ny+1; i++){
        for(int j=0; j<fmm.nx+1;j++){
            fmm.T_1D.at(a) = fmm.T.at(i).at(j);
            a++;
        }
    }

    fmm.export_vtu("result.vtu");

}