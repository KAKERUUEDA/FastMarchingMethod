#include "../lib/lib.h"
#include <sys/stat.h>
using namespace std;

int main()
{
    FMM fmm;

    std::string outputDir = "../../output";
    mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

    fmm.DefineGrid();
    fmm.FastMarchingMethod();
    fmm.updateT1D();

    std::string output = "../../output/result.vtu";
    fmm.export_vtu(output);
}