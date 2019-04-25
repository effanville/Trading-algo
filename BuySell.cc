#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <string>


#include <iostream>
#include <iomanip> 
#include <sstream>
#include <fstream>


#include "dayencoding.h"
#include "StockDecision.h"
#include "BuySell.h"
StPos::StPos(std::string name1, int n1,int Date1, double InitVal1)
{
    name=name1;
    n=n1;
    Date=Date1;
    InitVal=InitVal1;
};

std::string StPos::getname(void)
{
    return name;
};
int StPos::getn(void)
{
    return n;
}

int StPos::getDate(void)
{
    return Date;
};
double StPos::getInitVal(void)
{
    return InitVal;
};

double FindVal(std::vector<StPos> &Val,std::vector<Stock> &Info, std::string name,int Date)
{
    double sum=0;
    int i=0;
    for(i=0;i<Val.size();i++)
    {
        if(Val[i].getname()==name)
        {
            
            sum+=StockLatCl(Val[i].getname(),Info,Date)*Val[i].getn();
            
        }
        
    }
    return sum;
};

double Portfovalue(std::vector<StPos> &Val , std::vector<Stock> &Info,int t)
{
    double sum=0;
    int i=0;
    for(i=0;i<Val.size();i++)
    {
        sum+=StockLatCl(Val[i].getname(),Info,t)*Val[i].getn();
        
    }
    return sum;
};

/*int ReadPortfodata(double &A)
{
    std::ifstream Datainput;
    Datainput.open("Datainp.txt",std::ios::in);
    while(Datainput.good())
    {
        std::string line;
        getline(Datainput,line);
        std::string inp;
        inp=line.substr(0,4);
        if(inp=="frca")
        {
            std::string dat;
            dat=line;
            dat=line.erase(0,5);
            std::stringstream liner(dat);
            double num;
            while(liner >>num)
            {
                A=num;
            }
        }
    }
    return 1;
};

int WritePortfodata(double &A)
{
    std::ofstream WriteData;
    WriteData.open("Datainp.txt",std::ios::out);
    WriteData<< "frca " << A<<'\n';
    WriteData.close();
    return 1;
};




int ReadPortfo(std::vector<StPos> &Portfolio)
{
    std::ifstream CurPortfo;
    CurPortfo.open("Current.txt",std::ios::in);
    while(CurPortfo.good())
    {
        std::string line; //creates empty string
        getline(CurPortfo, line);
        std::string name;

        std::size_t found = line.find(" ");

        name = line.substr(0,found);

        line.erase(0,found+1);
        std::stringstream liner(line);
        double num;
        std::vector<double> info;
        while(liner>>num)
        {
            info.push_back(num);
        }
        if(info.size()==3)
        {
            int a=(int)info[0];
            int b=(int)info[1];
            StPos thisone(name,a,b,info[2]);
            Portfolio.push_back(thisone);
        }
    }
    CurPortfo.close();
    return 1;
};

int WritePortfo(std::vector<StPos> &Portfolio)
{
    std::ofstream CurPortfoOut;
    CurPortfoOut.open("Current.txt",std::ios::out);
    int i=0;
    for(i=0;i<Portfolio.size();i++)
    {
        CurPortfoOut<< Portfolio[i].getname()<< " " <<  Portfolio[i].getn()<< " "<< Portfolio[i].getDate()<< " "<< Portfolio[i].getInitVal()<< '\n';
    }
    CurPortfoOut.close();
    return 1;
};*/

int BuySellAlgo(double Init, double &freecash,std::vector<StPos> &Portfolio, std::vector<Stock> &Info, int date,double cost,double frac1,int dayahead,int &Scount,int &Bcount,double &AvheldU,double &AvheldD, int &LossNo, double &LossAmount,int &GainNo, double &GainAmount)
{
    
    //One cycles through the stocks listed for selling and buying and then, if enough cash is available to buy, one purchases the stock
    //Assume one buys as the highest value that day,
    //and assume one sells at the lowest value for that day.
    //One also needs enough available cash to exact all the buy sales
    
    //The following creates the total value of the portfolio at current date
    //This provides an 
    double Total = freecash+Portfovalue(Portfolio,Info,date);

    
    ///////////////////////////////////////////////////////////////////////////////////////////
    //This process reads in the lists of buy stocks and sell stocks from the decision algorithm
    std::vector<std::string> Buy;
    std::vector<std::string> Sell;
    std::ifstream SellList;
    std::string filename="Outputs//";
    filename.append("Sell");
    filename+=std::to_string(cost);
    filename.append("-");
    filename+=std::to_string(frac1);
    filename.append("-");
    filename+=std::to_string(dayahead);
    filename.append(".txt");
    
    SellList.open(filename.c_str(), std::ios::in);
    if (!SellList) 
    {//ensures that the file can indeed be opened.
        std::cout << "==============================================================="<<'\n';
        std::cout << " Unable to open file Sell.txt" << '\n';
        std::cout << "==============================================================="<<'\n';
    }
    while(SellList.good())
    {
        std::string line; //creates empty string
        getline(SellList, line); 
        if(line.empty()==false)
        {
            Sell.push_back(line);
        }
    }
    SellList.close();
    std::ifstream BuyList;
    
    filename="Outputs//";

    filename.append("Buy");
    filename+=std::to_string(cost);
    filename.append("-");
    filename+=std::to_string(frac1);
    filename.append("-");
    filename+=std::to_string(dayahead);
    filename.append(".txt");

    BuyList.open(filename.c_str(), std::ios::in);
    if (!BuyList) 
    {//ensures that the file can indeed be opened.
        std::cout << "==============================================================="<<'\n';
        std::cout << " Unable to open file Buy.txt" << '\n';
        std::cout << "==============================================================="<<'\n';
    }
    while(BuyList.good())
    {
        std::string line; //creates empty string
        getline(BuyList, line); 
        if(line.empty()==false)
        {
            Buy.push_back(line);
        }
    }
    BuyList.close();
    ///////////////////////////////////////////////////////////////////////////////////////////
    //This part executes the buy and sell algorithm itself
    int i=0;
    std::ofstream TransList;
    
    filename="Outputs//";
    filename.append("Transactions");
    filename+=std::to_string(cost);
    filename.append("-");
    filename+=std::to_string(frac1);
    filename.append("-");
    filename+=std::to_string(dayahead);
    filename.append(".txt");
    
    TransList.open(filename.c_str(), std::ios_base::app | std::ios_base::out);
    
    int k=-1;

    //std::cout<<"about to sell"<<'\n';
    for(i=0;i<Sell.size();i++)
    {
        int j=0;
        for(j=0;j<Portfolio.size();j++)
        {
            if(Portfolio[j].getname()==Sell[i])
            {
                //Sell Portfolio[k].n lots of shares of stock//
                //this price is modified by cost to take into account the inability to buy at the actual open price
                //this cost is either helpful or a hindrance based on a random number generator.

                /* generate secret number between 1 and 100: */
                int iS;
                iS = rand() % 100;
                int direction=0;
                if(iS>49)
                {
                    direction=1;
                    
                }
                if(iS<=50)
                {
                    direction=-1;
                    
                }
                
                //std::cout<<"SELL "<< iS<< " " << direction<<'\n';
                double sellprice = StockLatOp(Portfolio[j].getname(),Info,date)*(1+cost*((double)direction));
                //here look into future to set buy price. dangerous
                //double sellprice =StockLatOp(Portfolio[j].getname(),Info,date)-cost*(StockLatOp(Portfolio[j].getname(),Info,date)-StockLatLow(Portfolio[j].getname(),Info,date));
                
                TransList<< date <<" Sel "<< Sell[i]<< " " << Portfolio[j].getn()<<" "<< sellprice<<" "<<Portfolio[j].getn()*  sellprice<<'\n';

                
                freecash+=Portfolio[j].getn()*  sellprice;
                
                ////////////////////////////////////////////////////////////
                //The following are for calculating statistics about the algorithm only
                if(Portfolio[j].getn()*(sellprice- Portfolio[j].getInitVal())>=0)
                {
                    GainNo++;
                    GainAmount+=( sellprice-Portfolio[j].getInitVal())/(Portfolio[j].getInitVal());
                    AvheldU+=(double)Portfolio[j].getDate()-(double)date;
                }
                if(Portfolio[j].getn()*( sellprice-Portfolio[j].getInitVal())<0)
                {
                    LossNo++;
                    LossAmount+= (Portfolio[j].getInitVal() -sellprice )/(Portfolio[j].getInitVal());
                    AvheldD+=(double)Portfolio[j].getDate()-(double)date;
                }
                Scount++;
                /////////////////////////////////////////////////////////////////////
                
                Portfolio.erase(Portfolio.begin()+j);
                freecash-=6; //selling has a cost
            }
        }
        k=-1;

    }
    i=0;
    //std::cout<<"about to buy"<<'\n';
    for(i=0;i<Buy.size();i++)
    {

        if(freecash>600)
        {
            /* generate secret number between 1 and 100: */
            int iS;
            iS = rand() % 100;
            int direction=0;
            if(iS>49)
            {
                direction=1;
                    
            }
            if(iS<=50)
            {
                direction=-1;
                    
            }
            
            //std::cout<<"BUY "<< iS<< " " << direction<<'\n';
            double u = StockLatOp(Buy[i],Info,date)*(1+cost*((double)direction));
            
            //the following  is the a priori max we have to pay to buy, so is the amount we must have
            //to be able to buy this stock
            //i.e. we assume that stock will move a little in a bad direction for us
            //and so ensure we have enough to buy this
            double u2= StockLatOp(Buy[i],Info,date)*(1+cost);
           
            
            //assume cost one has to pay is somewhere between the open and high price
            //double u=(StockLatOp(Buy[i],Info,date)-cost*(StockLatOp(Buy[i],Info,date)-StockLatHigh(Buy[i],Info,date))); 
            int n=0;
            double a=0;
            //if(Info[i].VolTest(2,14,date)+Info[i].ADXTest(14,date)==2)
            //{
                
            //fractions in these two are up to decision
            a =std::min(std::min(Total*frac1,freecash),20000.0);
            //}
            //if(Info[i].VolTest(2,14,date)+Info[i].ADXTest(14,date)!=2)
            //{
            //    a =std::min(std::min(Total*frac1/2,freecash),20000.0);
            //}

            if(u!=0)
            {
                while(n*u<=a)
                {
                    n++;
                }
                
                n--;
                if(n>1)
                {
 
                    if(freecash>n*u2+6)
                    {

                        //Buy n shares of stock //;
                        StPos thatone(Buy[i],n,date,u);
                        Portfolio.push_back(thatone);
                        TransList<< date <<" Buy "<< Buy[i]<< " " << n<< " "<< u<<" "<<n*u<<'\n';
                        freecash-=n*u;
                        freecash -=6;
                        Bcount++;
                    }
                }
            }
            
        }
    }
    TransList.close();
    //std::cout<< "bought"<<'\n';
    return 1;
};
