#pragma once

#ifndef StockDecision_H
#define StockDecision_H

double Mn(std::vector<double> A);
double STD(std::vector<double> A );

class Stock
{
public:
    Stock(std::string name, std::vector<int> Times1, std::vector<double> I, std::vector<double> H, std::vector<double> L, std::vector<double> C, std::vector<double> V);
    std::vector<double> Trend(int n,std::vector<std::vector<double> > Spec,int lambda,int t);
    double PredVal(int n,std::vector<std::vector<double> >Spec, int lambda,int t);
    int TrendTest(int n,std::vector<std::vector<double> > Spec,int lambda,int t);
    
    bool VolTest(int n,int m,int t);
    std::vector<double> DMplus(int n,int t);
    std::vector<double> DMminus(int n,int t);
    std::vector<double> TR(int n,int t);
    std::vector<double> DIplus(int n,int t);
    std::vector<double> DIminus(int n,int t);
    double ADX(int n,int t);
    bool ADXTest(int n,int t);
    std::vector<double> Stochastic(int n,int t);
    int StochTest(int n,int t);
    int MACTest(int n, int m,int t);
    void PriceUpdate(double I2, double H2, double L2,double C2, double V2);
    int CandTest(int n,std::vector<std::vector<double> >Spec,int lambda,int t);
    int BuySellDec(int t,std::vector<std::vector<double> >Spec,int nTr,int lambda, int nV,int mV,int nA, int nS1,int nS2,int nMAC1,int mMAC1, int nMAC2,int mMAC2,int dayahead);
    std::string nameOut(void);
    std::vector<double> COut(void);
    std::vector<double> IOut(void);
    std::vector<double> HOut(void);
    std::vector<double> LOut(void);
    std::vector<double> VOut(void);
    int HistLength(void);
    std::vector<int> Times;
    std::vector<double> BetaVals; //stores the predictions of beta for each stock
    
    double LassoCV(int m,int nTr,std::vector<std::vector < double> > &Spec,int starttime,int lambda2, int nV,int mV,int nA,int nS1,int nS2,int nMAC1,int mMAC1,int nMAC2,int mMAC2,int dayahead);
    
    double errorF(double lambda,int m,std::vector<int> Fk,int nTr,std::vector<std::vector < double> > &Spec,int starttime,int lambda2, int nV,int mV,int nA,int nS1,int nS2,int nMAC1,int mMAC1,int nMAC2,int mMAC2,int dayahead);//produces error for predictor not including Fk and for lambda 
    
    std::vector<double> LPs(double lambda,int m,int nTr, std::vector<std::vector < double> > &Spec,int starttime,int lambda2, int nV,int mV,int nA,int nS1,int nS2,int nMAC1,int mMAC1,int nMAC2,int mMAC2,int dayahead);//creates the prediction for each stock
    
    double LPred(int t, int nTr,std::vector<std::vector < double> > &Spec,int lambda2, int nV,int mV,int nA,int nS1,int nS2,int nMAC1,int mMAC1,int nMAC2,int mMAC2,int dayahead);//gives a value as output for future 5 days

private:
    std::string name;
    std::vector<double> I;//opening price
    std::vector<double> H;//daily high price
    std::vector<double> L;//daily low price
    std::vector<double> C;//closing price
    std::vector<double> V; //volume of that day
};


double MSR(std::vector<double> A, std::vector<double> B);

float norm_feature(int j, std::vector<std::vector<double> >arr, int n);
float approx(int dim,int d_ignore,std::vector<double> weights,std::vector<std::vector<double> >arr,int 
i);
float rho_j(std::vector<std::vector<double> >arr,int n,int j,int num_dim,std::vector<double> weights);
float intercept(std::vector<double> arr1,std::vector<std::vector<double> >arr,int num_dim,int num_obs);

std::vector<double> Lasso(const int num_dim,int num_obs,std::vector<std::vector<double> > & data,std::vector<double> &y,double lambda=0.1);

//calculates moving average of vector v starting at index t for the following n elements
double MA(std::vector<double>v,int n,int t);

double EMA(std::vector<double>v,double alpha,int n,int t);

std::vector<double> LassoFit(std::vector<Stock> &Allstock,double lambda, int n,int m);

double StockLatHigh(std::string name, std::vector<Stock> &St,int time);

double StockLatLow(std::string name, std::vector<Stock> &St,int time);

double StockLatOp(std::string name, std::vector<Stock> &St,int time);

double StockLatCl(std::string name, std::vector<Stock> &St,int time);

//////////////////////////////////////////////////////////

//for outputting and inputting data

int StockInput(std::vector<Stock> &Input);
int NStockInput(std::vector<Stock> &Input);
int WriteStockList(std::vector<Stock> Input);

//produces a list of all possible stocks with their buy sell data
void StockOutput(std::vector<Stock> &Universe,int nTr,std::vector<std::vector<double> > &Spec,int lambda, int nV,int mV,int nA,int nS1,int nS2,int nMAC1,int mMAC1,int nMAC2,int mMAC2,int t,double FC,double Cost,double frac,int dayahead);

//produces a file with a collection of stocks to buy
void StockOutputBuy(std::vector<Stock> &Universe,int nTr,std::vector<std::vector<double> > &Spec,int lambda, int nV,int mV,int nA,int nS1, int nS2,int nMAC1,int mMAC1,int nMAC2,int mMAC2,int t,double FC,double Cost,double frac,int dayahead);
//produces a file with a list of stocks to sell
void StockOutputSell(std::vector<Stock> &Universe,int nTr,std::vector<std::vector<double> > &Spec,int lambda, int nV,int mV,int nA,int nS1,int nS2,int nMAC1,int mMAC1,int nMAC2,int mMAC2,int t,double FC,double Cost,double frac,int dayahead);

#endif
