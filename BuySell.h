#pragma once

#ifndef BuySell_H
#define BuySell_H
//the following is a class to describe how many shares one owns of a specific stock


class StPos
{
public:
    StPos(std::string name, int n, int Date, double InitVal);
    std::string getname(void);
    int getn(void);
    int getDate(void);
    double getInitVal(void);
private:
    std::string name;
    int n;//number of shares bought
    int Date;//integer representing the date stock was bought on
    double InitVal;//value of one share at time of buying
};

double FindVal(std::vector<StPos> &Val ,std::vector<Stock> &Info, std::string name,int Date);

int ReadPortfodata(double &A);

int WritePortfodata(double &A);

double Portfovalue(std::vector<StPos> &Val , std::vector<Stock> &Info,int t);

int ReadPortfo(std::vector<StPos> &Portfolio);

int WritePortfo(std::vector<StPos> &Portfolio);

int BuySellAlgo(double Init,double &freecash,std::vector<StPos> &Portfolio, std::vector<Stock> &info, int date,double cost,double frac,int dayahead, int &Scount,int &Bcount,double &AvheldU,double &AvheldD, int &LossNo, double &LossAmount,int &GainNo, double &GainAmount);

#endif
