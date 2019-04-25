#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <string>
#include <time.h>



#include <iostream>
#include <iomanip> 
#include <sstream>
#include <fstream>
#include "quicksort.h"
#include "dayencoding.h"

#include "StockDecision.h"

#include "BuySell.h"
///////////////////////////////////////////////////////
//function returns compound annual rate.
//inputs should be added in order of time1, price1, time2, price2 where time1>time2 
double YrAvRet(std::vector<double> A)
{
    double sum=0;
    int i =0;
    for(i=1;i<A.size();i++)
    {
        sum+=(A[i]-A[i-1])/A[i];
    }
    return sum/(A.size()-1);
};

double YrAvVar(std::vector<double> A)
{
    int i=0;
    double sum=0;
    double m=YrAvRet(A);
    for(i=1;i<A.size();i++)
    {
        sum+=pow(((A[i]-A[i-1])/A[i-1])-m,2);
    }
    return sum/(A.size()-1);
};



//////////////////////////////////////////////////
double CARab (int a, double x, int b, double y){
    if((y==0)&& (x==0))
    {
        return 0;
    }
    if(y==0)
    {
        return 0;
    }
    int c =a-b;

    double d= 365/(double) c;
    return (pow(x/y,d)-1);  
}; 


///////////////////////////////////////////////
//This is used for inverting matrix
static int LUPdecompose(int size, std::vector< std::vector<double> > &A, std::vector<int> &P)
   {
    int i, j, k, kd = 0, T;
    double p, t;

 /* Finding the pivot of the LUP decomposition. */
    for(i=0; i<size; i++) 
    {
        P[i]=i;
    } //Initializing.
    
    for(k=0; k<size-1; k++)
       {
        p = 0;
        t = A[k][k];
        
        if(t!=0)
        {
            p=t;
        }

        if(p == 0)
           {
            printf("\nLUPdecompose(): ERROR: A singular matrix is supplied.\n"\
                   "\tRefusing to proceed any further.\n");
            return -1;
           }

        for(i=k+1; i<size; i++) //Performing substraction to decompose A as LU.
            {
             A[i][k] = A[i][k]/A[k][k];
             for(j=k+1; j<size; j++) 
             {
                 A[i][j] -= A[i][k]*A[k][j];
             }
            }
        } //Now, 'A' contains the L (without the diagonal elements, which are all 1)
          //and the U.
    return 0;
   }




/* This function calculates the inverse of the LUP decomposed matrix 'LU' and pivoting
 * information stored in 'P'. The inverse is returned through the matrix 'LU' itselt.
 * 'B', X', and 'Y' are used as temporary spaces. */
int LUPinverse(int size, std::vector< int> P, std::vector< std::vector<double> > &LU )
   {
       std::vector< std::vector<double> > B;
       B.resize(size);
       std::vector< double> X;
       X.resize(size);
       std::vector<double> Y;
       Y.resize(size);
    int i, j, n, m;
    for(i=0;i<size;i++)
    {
        B[i].resize(size);
    }
    double t;
 /* Solving LUX = Pe, in order to calculate the inverse of 'A'. Here, 'e' is a column
  * vector of the identity matrix of size 'size-1'. Solving for all 'e'. */
    for(i=0; i<size; i++)
     {
    //Storing elements of the i-th column of the identity matrix in i-th row of 'B'.
      for(j = 0; j<size; j++)
      {
          B[i][j] = 0;
      }
      B[i][i] = 1;
     }
   //Solving Ly = Pb.
   for(i=0; i<size; i++)
     {
     for(n=0; n<size; n++)
       {
        t = 0;
        for(m=0; m<n; m++) 
        {
            t += LU[n][m]*Y[m];
        }
        Y[n] = B[P[n]][i]-t;
       }
   //Solving Ux = y.
     for(n=size-1; n>=0; n--)
       {
        t = 0;
        for(m = n+1; m < size; m++)
        {
            t += LU[n][m]*X[m];
        }
        X[n] = (Y[n]-t)/LU[n][n];

       }//Now, X contains the solution.

      for(j = 0; j<size; j++)
      {
          B[j][i] = X[j]; //Copying 'X' into the same row of 'B'.
     } //Now, 'B' the transpose of the inverse of 'A'.
     }
 /* Copying transpose of 'B' into 'LU', which would the inverse of 'A'. */
    for(i=0; i<size; i++) for(j=0; j<size; j++) LU[i][j] = B[i][j];

    return 1;
   }

//////////////////////////////////////////////////////////

//Setting default parameters

double FC=20000.0;
int end=-1;
double Init=FC;
double cost=0.0;
double frac=0.1;
int dayahead=2;


///////////////////////////////////////////////////////////////////

int main(int argc,char *argv[] )
{

    //The following takes user input parameters and adds to the program
    //If none specified program takes the defaults as above
    if( argc ==4 ) 
    {
        cost = atof(argv[1]);  
        frac = atof(argv[2]);  
        dayahead =atoi(argv[3]);
    }
    
    int i=0;
    srand(time(NULL));
    std::cout << " Cost is " << cost<<'\n';
    std::cout << " Frac is " <<frac<<'\n';
    std::cout << " Prediction Day is " <<dayahead<<'\n';

    /////////////////////////////////////////////////
    //Closk is used for timing purposes only
    clock_t tStart = clock();

    
    ////////////////////////////////////////////////////////////////////
    //The following is the algorithm that reads in data on the stock universe into stock classes for use throughout the remainder of the algorithm

    std::vector<Stock> AllStocks;
    
    int m= NStockInput(AllStocks);

    //It is necessary for the data to be displayed in order of newest value to oldest value
    if(m==-1)
    {
        std::cout <<  "Could not populate database"<<'\n';
        return 1;
    }
    std::cout << "Created Stock Database"<<'\n';
    std::cout << "have " << m<< " stocks"<<'\n';
    //The "time" is the index of the data values
    //the time thus gets smaller and smaller until the present, which is time=0
    
    int start=20;
    int nm=0;
    int time=AllStocks[0].HistLength();
    int index=0;
    for(nm=0;nm<AllStocks.size();nm++)
    {
        if(AllStocks[nm].HistLength()<time)
        {    
            time=AllStocks[nm].HistLength();
            index =nm;
        }
    }
    time =time-start;
    
    int starttime = time;
    std::cout << "lowest number with "<<AllStocks[index].nameOut()<<" "<<AllStocks[index].HistLength()<<'\n';    

    //this needs to be done once
    int n=10;
    int lambda=5000;
    int nV=2;
    int mV=14;
    int nA=14;
    int nS1=14;
    int nS2=28;
    int nS3=7;
    int nMAC1=5;
    int mMAC1=20;
    int nMAC2=20;
    int mMAC2=50;
    int nMAC3=50;
    int mMac3=200;
    
    
    ///////////////////////////////////////////////////////////
    //calculation for hodrick prescot filter
    std::vector<std::vector<double> > Spec(n, std::vector<double>(n,0.0));  
    std::vector<std::vector<double> > K(n, std::vector<double>(n,0.0)); 
    if(n>6)
    {
        K[0][0]=1;
        K[1][1]=5;
        K[1][0]=-2;
        K[0][1]=-2;
        K[0][2]=1;
        K[1][2]=-4;
        K[1][3]=1;

        for(i=2;i<K.size()-2;i++)
        {
            K[i][i]=6;
            K[i][i-1]=-4;
            K[i][i-2]=1;
            K[i][i+1]=-4;
            K[i][i+2]=1;
            
        }
        K[K.size()-1][K.size()-1]=1;
        K[K.size()-2][K.size()-2]=5;
        K[K.size()-2][K.size()-1]=-2;
        K[K.size()-1][K.size()-2]=-2;
        K[K.size()-1][K.size()-3]=1;
        K[K.size()-3][K.size()-2]=-4;
        K[K.size()-2][K.size()-3]=-4;
        K[K.size()-2][K.size()-4]=1;
    }

    int j=0;
    for(i=0;i<Spec.size();i++)
    {
        for(j=0;j<Spec.size();j++)
        {
            if(i==j)
            {
                Spec[i][j]=1+lambda*K[i][j];
            }
            if(i!=j)
            {
               Spec[i][j]=lambda*K[i][j]; 
            }
        }
    }
        
   std::vector<int> P(Spec.size(),1.0);

    LUPdecompose(Spec.size(),Spec, P); 

    LUPinverse(Spec.size(), P, Spec);
    
    std::vector<int>().swap(P);


    
    
    double benchval=0;
    double oldbenchval=0;
    double lastportfoval=0;

    
    /////////////////////////////////////////////////////////////
    //here we calculate all the prediction values for each stock
    //inputs are 
    //double lambda,int m, int nTr,std::vector<std::vector < double> > Spec,int starttime,int lambda2, int nV,int mV,int nA,int nS1,int nS2,int nMAC1,int mMAC1,int nMAC2,int mMAC2,int dayahead
    int periods=100;
    for(i=0;i<AllStocks.size();i++)
    {
        double lammmy=0.1;
        lammmy=AllStocks[i].LassoCV(periods, n,Spec, starttime, lambda, nV, mV, nA, nS1, nS2, nMAC1, mMAC1, nMAC2, mMAC2,dayahead);
        
        AllStocks[i].BetaVals=AllStocks[i].LPs(lammmy,periods, n,Spec, starttime, lambda, nV, mV, nA, nS1, nS2, nMAC1, mMAC1, nMAC2, mMAC2,dayahead);
    }
    //set time to the end of the fitting period
    time =time- periods*20;
    starttime=time;
    
                /////////////////////////////////
    //for the purposes of benchmark portfolio
    //This is a flat portfolio of equal weighting on all stocks in portfolio
    double ratioval = Init/((double)m); 
    std::vector<double> ratios(m,0.0);
    int k=0;
    for(k=0;k<AllStocks.size();k++)
    {
        ratios[k]=ratioval/StockLatOp(AllStocks[k].nameOut(),AllStocks,starttime);
    }
    std::cout<< "produced betas"<<'\n';
    
    //while(frac<=0.5)
    //{    
        std::vector<double> yearvals;
        std::vector<StPos> Portfo;
            time=starttime;
            FC=20000.0;
            
    ////////////////////////////////////////////////////////////
    //These are parameters used solely for information on the test run
    
    int BuyCount=0;
    int SellCount=0;
    double LossA=0;
    int LossN=0;
    double GainA=0;
    int GainN=0;
    double AvHeldLenU=0;
    double AvHeldLenD=0;
    double AvRelInv=0;
    int AvHoldingNo=0;
    
    /////////////////////////////////////////////////////////////
    //initialises output stream for portfolio data
    std::ofstream PortfoCur;
    std::string filename="Outputs//";
    filename.append("CurVal");
    filename+=std::to_string(cost);
    filename.append("-");
    filename+=std::to_string(frac);
    filename.append("-");
    filename+=std::to_string(dayahead);
    filename.append(".txt");
    
    PortfoCur.open(filename.c_str(), std::ios_base::app | std::ios_base::out);
    PortfoCur << "Simulation with "<<m<< " stocks starting on "<<starttime<<  '\n';
    PortfoCur << "Fitting started on "<<AllStocks[index].HistLength()<<'\n';
    PortfoCur << " Cost is " << cost<<'\n';
    PortfoCur << " Frac is " <<frac<<'\n';
    PortfoCur << " Prediction Day is " <<dayahead<<'\n';
    
    PortfoCur<< "date    ";
    for(i=0;i<AllStocks.size();i++)
        {
            PortfoCur<< AllStocks[i].nameOut()<<" ";
        }
    PortfoCur <<"FC ";
    PortfoCur<<'\n';
    /////////////////////////////////////////////////////
    
    while(time>end)
    {
        //The following is the Decision algorithm that decides on what to buy and what to sell at the current times
        
        //This work is implicit. It prints to files Buy.txt and Sell.txt the decisions
        //These are all that are needed for the execution algorithm
        //The decisions are contained in the function BuySellDec within these
        //The parameters are standard and can be found in "trading strategy.pdf"
        //under the relevant definitions for the various
        //The order is
        //int nTr,int lambda, int nV,int mV,int nA, int nS,int nMAC1,int mMac1, int nMAC2,int mMAC2,int t
        //the integer t denotes the time, so currently 0 means current time

        StockOutput(AllStocks,n,Spec,lambda, nV, mV, nA, nS1, nS2, nMAC1, mMAC1, nMAC2, mMAC2, time,Init,cost,frac,dayahead);
            
        ///////////////////////////////////////////////////////////////////////////////////////////
        //The following is the execution algorithm. 
        //Requires files Buy.txt and Sell.txt containing lists of the names of stocks to buy and sell

        
        //ReadPortfo(Portfo); //necessary for live trading where run program for short time each day
        //Executes the buy and sell orders in Buy.txt and Sell.txt, updating the portfolio as it proceeds
        //Also writes out the transactions performed into file Transactions.txt
        
        // includes argument time
        //note one can set this to today for execution, or iterate through dates for testing algorithm purposes.

        BuySellAlgo(Init,FC,Portfo, AllStocks, time, cost,frac,dayahead, SellCount, BuyCount, AvHeldLenU,AvHeldLenD,LossN, LossA,GainN, GainA);

        //Writes the current portfolio to file Current.txt
        //WritePortfo(Portfo);
        
        /////////////////////////////////////////////////////////////////////////////
        //The remaining part of code is just for creating statistics of the simulation
        
            //calculates average relative investment of cash
        AvRelInv+=Portfovalue(Portfo,AllStocks,time)/(FC+Portfovalue(Portfo,AllStocks,time));
        AvHoldingNo+=Portfo.size();
        
        /////////////////////////////////////
        
        PortfoCur<<time << " "<<Dayinv(AllStocks[index].Times[time])[0]<<" "<<Dayinv(AllStocks[index].Times[time])[1] << " " <<Dayinv(AllStocks[index].Times[time])[2]<<" ";
        
        oldbenchval=benchval;
        
        benchval=0;
        
        for(i=0;i<AllStocks.size();i++)
        {
            PortfoCur<< FindVal(Portfo, AllStocks, AllStocks[i].nameOut(), time)<<" ";
            benchval+=StockLatOp(AllStocks[i].nameOut(),AllStocks,time)*ratios[i];
        }
        
        PortfoCur << FC<<" "<<FC+Portfovalue(Portfo,AllStocks,time)<<" " << benchval<<" " << AvRelInv/((double)starttime-(double)time)<<" " << AvHoldingNo/((double)starttime-(double)time)<<'\n';
        if((time-starttime)%250==0)
        {
            yearvals.push_back(FC+Portfovalue(Portfo,AllStocks,time));
        }
            
        if(time % 100==0)
        {
            std::cout <<time <<" "<<Dayinv(AllStocks[index].Times[time])[0] <<" "<< Dayinv(AllStocks[index].Times[time])[1] << " "  <<Dayinv(AllStocks[index].Times[time])[2]<<" "<<FC << " " << FC+Portfovalue(Portfo,AllStocks,time) << " " << benchval<<" " << AvRelInv/((double)starttime-(double)time)<<" "<< AvHoldingNo/((double)starttime-(double)time)<<'\n';
        }
        
        lastportfoval=FC+Portfovalue(Portfo,AllStocks,time);
        
        //////////////////////////////////////////////////////////////////////////////
        //replaces time with the next value to be used
        time--;
    }

    
    //Prints out information about the run of the test to the relevant file and also to the command line
    PortfoCur<< '\n';
    PortfoCur << "BuyCount SellCount  AvHeldLength"<<'\n';    
    PortfoCur << BuyCount<<" " <<SellCount<< " " <<(AvHeldLenU+AvHeldLenD)/BuyCount<<'\n';
    PortfoCur<< "NoLosses AvLossAmount AvHeldDLength Fraction"<<'\n';    
    PortfoCur << LossN<<" " <<LossA/LossN<< " "<< AvHeldLenD/(double)LossN <<" " <<(double)LossN/((double)LossN+(double)GainN)<<'\n';
    PortfoCur << "NoGains AvGainAmount AvHeldULength Fraction"<<'\n';    
    PortfoCur << GainN<<" " <<GainA/GainN<< " " << AvHeldLenU/(double)GainN<<" "<<(double)GainN/((double)LossN+(double)GainN)<<'\n';
    PortfoCur << "CAR " << CARab(AllStocks[index].Times[time+1],FC+Portfovalue(Portfo,AllStocks,time+1),AllStocks[index].Times[starttime],Init)<<'\n';
    PortfoCur << YrAvRet(yearvals)<< " "<< YrAvVar(yearvals)<< " Sharpe " << YrAvRet(yearvals)/(pow( YrAvVar(yearvals),0.5))<<'\n';
    
    PortfoCur.close();

    std::cout << "BuyCount SellCount AvHeldLength"<<'\n';    
    std::cout << BuyCount<<" " <<SellCount<< " " <<(AvHeldLenU+AvHeldLenD)/BuyCount<<'\n';
    std::cout << "NoLosses AvLossAmount AvHeldDLength Fraction"<<'\n';    
    std::cout << LossN<<" " <<LossA/LossN<< " "<< AvHeldLenD/(double)LossN <<" " <<(double)LossN/((double)LossN+(double)GainN)<<'\n';
    std::cout << "NoGains AvGainAmount AvHeldULength Fraction"<<'\n';    
    std::cout << GainN<<" " <<GainA/GainN<< " " << AvHeldLenU/(double)GainN<<" "<<(double)GainN/((double)LossN+(double)GainN)<<'\n';
    std::cout << "CAR " << CARab(AllStocks[index].Times[time+1],FC+Portfovalue(Portfo,AllStocks,time+1), AllStocks[index].Times[starttime] ,Init)<<'\n';
    std::cout << YrAvRet(yearvals)<< " "<< YrAvVar(yearvals)<< " Sharpe " << YrAvRet(yearvals)/(pow( YrAvVar(yearvals),0.5))<<'\n';
    
        //frac=frac+0.1;
    //}
    
    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    
    return 0;
}
