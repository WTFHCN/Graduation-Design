#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <assert.h>
double calcError(std::vector<double> &reality, std::vector<double> &ans)
{
    assert(reality.size() == ans.size());
    double res = 0;
    for (int i = 0; i < ans.size(); i++)
    {
        res += (reality[i] - ans[i]) * (reality[i] - ans[i]);
    }
    return sqrt(double(res / reality.size()));
}
void wirteCsv(std::string &fileName, std::vector<std::string> &name, std::vector<std::vector<double>> &data)
{
    std::ofstream outFile;
    assert(name.size() == data.size());
    outFile.open(fileName, std::ios::out);
    for (int i = 1; i < name.size(); i++)
    {
        outFile << name[i];
        // for (int j = 0; j < data[i].size(); j++)
        // {
        outFile << std::setprecision(6);
        outFile << ',' << calcError(data[0], data[i]);
        // }
        outFile << "\n";
    }
}