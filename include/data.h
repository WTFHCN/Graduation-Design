#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <assert.h>
#include "util.h"
void WirteCsv(std::string &fileName, std::vector<std::string> &name, std::vector<std::vector<double>> &error)
{
    std::ofstream outFile;
    // assert(name.size() == data.size());
    outFile.open(fileName, std::ios::out);
    for (int i = 3; i < error.size(); i++)
    {
        // std::std::cout << name[i] << '\n';
        outFile << name[i];
        for (int j = 8; j < error[i].size(); j++)
        {
            outFile << std::setprecision(6);
            outFile << ',' << error[i][j];
        }
        outFile << "\n";
    }
    outFile.close();
}

struct RetailData
{
    int StockCode;
    std::string Description;
    int Quantity;
    double UnitPrice;
    std::string Country;
};
std::vector<RetailData> ReadRetailCsv(std::string fileName = "Retail.csv")
{
    std::vector<RetailData> CensusList;
    std::ifstream fin(fileName);

    std::string line;
    getline(fin, line);
    while (getline(fin, line))
    {

        std::istringstream sin(line);
        // std::cout << line << '\n';
        std::string info;
        std::vector<std::string> infoList;

        while (getline(sin, info, ','))
        {
            infoList.push_back(info);
        }

        RetailData dataBase;
        dataBase.StockCode = atoi(infoList[2].c_str());
        dataBase.Description = infoList[2];

        dataBase.Quantity = atoi(infoList[3].c_str());
        dataBase.UnitPrice = atof(infoList[5].c_str());
        dataBase.Country = infoList[7].c_str();
        
        CensusList.push_back(dataBase);
    }
    return CensusList;
}
struct CensusData
{
    int age;
    int zipcode;
    int income;
    int vote;
    std::string sex, race, UNION;
};
std::vector<CensusData> ReadCensusCsv(std::string fileName = "census.csv")
{
    std::vector<CensusData> CensusList;
    std::ifstream fin(fileName);

    std::string line;
    while (getline(fin, line))
    {
        std::istringstream sin(line);
        std::string info;
        std::vector<std::string> infoList;

        while (getline(sin, info, ','))
        {
            infoList.push_back(info);
        }

        CensusData dataBase;
        dataBase.age = atoi(infoList[0].c_str());
        dataBase.sex = infoList[1];
        dataBase.race = infoList[2];
        dataBase.UNION = infoList[3];
        dataBase.zipcode = atoi(infoList[4].c_str());
        dataBase.income = atoi(infoList[5].c_str());
        dataBase.vote = atoi(infoList[6].c_str());

        CensusList.push_back(dataBase);
    }
    return CensusList;
}
void showCensus(std::vector<CensusData> &CensusList)
{
    for (int j = 0; j < CensusList.size(); j++)
    {
        auto dataBase = CensusList[j];
        std::cout << "database " << j << " age: " << dataBase.age << '\n';
        std::cout << "database " << j << " sex: " << dataBase.sex << '\n';
        std::cout << "database " << j << " race: " << dataBase.race << '\n';
        std::cout << "database " << j << " UNION: " << dataBase.UNION << '\n';
        std::cout << "database " << j << " zipcode: " << dataBase.zipcode << '\n';
        std::cout << "database " << j << " income: " << dataBase.income << '\n';
        std::cout << "database " << j << " vote: " << dataBase.vote << '\n';
        std::cout << '\n';
    }
}