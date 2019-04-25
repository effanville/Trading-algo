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

#include "quicksort.h"

extern Day Reference;
double Mn(std::vector<double> A)
{
    double sum=0;
    int i=0;
    for(i=0;i<A.size();i++)
    {
        sum+=A[i];
    }
    sum=sum/A.size();
    return sum;
}

double STD(std::vector<double> A )
{
    double sum=0;
    double m=Mn(A);
    int i =0;
    for(i=0;i<A.size();i++)
    {
        sum+=(A[i]-m)*(A[i]-m);
    }
    sum=sum/(A.size()-1);
    sum=pow(sum,0.5);
    return sum;
}
/////////////////////////////////////////
//functions for lasso regression
double MSR(std::vector<double> A, std::vector<double> B)
{
    if(A.size()!=B.size())
    {
        printf("Vectors Incompatible");
        printf("\n");
        return -1;
    }
    
    double sum=0;
    int i =0;
    for(i=0;i<A.size();i++)
    {
        sum+=pow(A[i]-B[i],2);
    }
    return sum;
}

/*Computes the normalization factor*/
float norm_feature(int j, std::vector<std::vector<double> >arr, int n)
{

    float sum = 0.0;
    int i;
    for(i=0;i<n;i++)
    {

        sum = sum + pow(arr[i][j],2);
    }

    return sum;
}


/*Computes the partial sum*/
float approx(int dim,int d_ignore,std::vector<double> weights,std::vector<std::vector<double> >arr,int 
i)
{

    int flag = 1;
    
    if(d_ignore == -1)
        flag = 0;
    
    int j;
    
    float sum = 0.0;
    
    for(j=0;j<dim;j++)
    {
        if(j != d_ignore)
            sum = sum + weights[j]*arr[i][j];
        else
            continue;
    }

    return sum;

}

/* Computes rho-j */
float rho_j(std::vector<std::vector<double> >arr,int n,int j,int num_dim,std::vector<double> weights)
{

    float sum = 0.0;
    int i;
    float partial_sum ;
    for(i=0;i<n;i++)
    {

        partial_sum = approx(num_dim,j,weights,arr,i);
        
        sum = sum + arr[i][j]*(arr[i][num_dim]-partial_sum);


    }

    return sum;
}

float intercept(std::vector<double> arr1,std::vector<std::vector<double> >arr,int num_dim,int num_obs) {

    int i;
    float sum =0.0;
    for (i = 0; i < num_obs; i++) 
    {
        sum = sum + pow((arr[i][num_dim]) -     approx(num_dim, -1, arr1, arr, i), 1);
    }

    return sum;

}






std::vector<double> Lasso(const int num_dim,int num_obs,std::vector<std::vector<double> > &data,std::vector<double> &y,double lambda)
{
    //create number of parameters and observations
    const int num_dim2 =num_dim+1;
    if(data.size()!=num_obs)
    {
        std::cout<< "not enough observations"<<'\n';
        std::vector<double> out;
        return out;
        
    }
    if(data[0].size()!=num_dim)
    {
        std::cout<< "not enough variables"<<'\n';
        std::vector<double> out;
        return out;
        
    }
    int i=0,j=0;
    //need to add the observation to the end of the data matrix for this to work
    for(i=0;i<num_obs;i++)
    {
        data[i].push_back(y[i]);
    }
    
    float a = 1.0;

    std::vector<double> weights(num_dim+1); //weights[end] contains the intercept
    std::vector<double> OldWeights(num_dim+1);


    //Initializing the weight vector
    for (i = 0; i < num_dim; i++)
        weights[i]=0.1;
        //weights[i] = ((float) rand() / (float) (RAND_MAX)) * a;

    int iter = 500;
    int t = 0;
    int r, l;
    std::vector<double> rho(num_dim);
    
    for (i = 0; i < num_dim; i++) 
    {
        rho[i] = rho_j(data, num_obs, r,num_dim, weights);
    }

    // Intercept initialization
    //weights[num_dim] = intercept(weights,data,num_dim, num_obs);
    weights[num_dim] = 0;

    while (t < iter) 
    {
        //OldWeights=weights;
        for (r = 0; r < num_dim; r++) 
        {
            rho[r] = rho_j(data, num_obs, r,num_dim, weights);
            //printf("rho %d:%f ",r,rho[r]);
            if (rho[r] < -lambda / 2)
            {
                weights[r] = (rho[r] + lambda / 2) / norm_feature(r, data, num_obs);
            }
            else if (rho[r] > lambda / 2)
            {
                weights[r] = (rho[r] - lambda / 2) / norm_feature(r, data, num_obs);
            }else
            {
                weights[r] = 0;
            }
            
            weights[num_dim] = 0;
            //intercept(weights, data, num_dim, num_obs);


        }
        if(MSR(OldWeights,weights) < 5e-4)
        {
            break;
        }

        t++;
    }

    
 return weights;

}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

double Stock::LassoCV(int m,int nTr,std::vector<std::vector < double> > &Spec,int starttime,int lambda2, int nV,int mV,int nA,int nS1,int nS2,int nMAC1,int mMAC1,int nMAC2,int mMAC2,int dayahead)
{
    
    int i =0;
    int n=10;
    int k=0;
    int K=5;
    int tot=m/K;
    
    std::vector<double> Errlam(n-1,0.0);
    std::vector<std::vector<int> > Fk(K, std::vector<int>(tot,0));
    for(i=0;i<K;i++)
    {
        for(k=0;k<tot;k++)
        {

            Fk[i][k]=tot*i+k;
        }
    }
    
    for(i=1;i<n;i++)
    {
        for(k=0;k<Fk.size();k++)
        {
            double lambda =(double)i/10.0;
            Errlam[i-1] += errorF(lambda, m,Fk[k],nTr, Spec, starttime, lambda2,  nV, mV, nA, nS1, nS2, nMAC1, mMAC1, nMAC2, mMAC2,dayahead);
        }
        
    }
    int maxElementIndex = std::max_element(Errlam.begin(),Errlam.end()) - Errlam.begin();
    return (maxElementIndex+1)/10.0;
}

//error here is a simple mean of difference

double Stock::errorF(double lambda,int m,std::vector<int> Fk,int nTr,std::vector<std::vector < double> > &Spec,int starttime,int lambda2, int nV,int mV,int nA,int nS1,int nS2,int nMAC1,int mMAC1,int nMAC2,int mMAC2,int dayahead)
{
    std::vector<double> y;
    std::vector<std::vector<double> >X;
    int t =starttime;
    int i=0;
    int j=0;
    for(j=0;j<m;j++)//cycles through all time periods
    {
        if(std::find(Fk.begin(), Fk.end(), j)==Fk.end())
        {
            t=starttime-j*20;
            std::vector<double> v;
            
            X.push_back(v);
            X.back().push_back(H[t]/I[t]);
            X.back().push_back(L[t]/I[t]);
            X.back().push_back(C[t]/I[t]);
            X.back().push_back(PredVal(nTr,Spec,lambda2,t));
            X.back().push_back(ADX(nA,t)/100.0);
            
            std::vector<double> Ied(mMAC2,0.0);

            for(i=0;i<mMAC2;i++)
            {
                Ied[i]=I[t+1+i]/I[t+1];
            }

            std::vector<double> Ved(mV,0.0);
            for(i=0;i<mV;i++)
            {
                Ved[i]=V[t+1+i]/V[t+1];
            }
            
            X.back().push_back(MA(Ied,nMAC1,0));
            X.back().push_back(MA(Ied,nMAC2,0));
            X.back().push_back(MA(Ied,mMAC2,0));

            X.back().push_back(MA(Stochastic(nS1,t),nS1,0)/100.0);
            X.back().push_back(MA(Stochastic(nS2,t),nS2,0)/100.0);
                    
            X.back().push_back(MA(Ved,nV,0));
            X.back().push_back(MA(Ved,mV,0));
                    
            y.push_back(I[t-dayahead]/I[t]);
        }
    }
    const int num_dim =X[0].size();
    int num_obs = X.size();
    std::vector<double> ouptup=Lasso(num_dim,num_obs,X, y,lambda);

    i=0;
    std::vector<double> preds;
    std::vector<double> actual;
    while(i<Fk.size())
    {
        int nt=starttime-Fk[i]*20;
        std::vector<double> Ied(mMAC2,0.0);

        for(j=0;j<mMAC2;j++)
        {
            Ied[j]=I[nt+1+j]/I[nt+1];
        }

        std::vector<double> Ved(mV,0.0);
        for(j=0;j<mV;j++)
        {
            Ved[j]=V[nt+1+j]/V[nt+1];
        }
        preds.push_back(ouptup[0]*H[nt+1]/I[nt+1]+ouptup[1]*L[nt+1]/I[nt+1]+ouptup[2]*C[nt+1]/I[nt+1]+ouptup[3]*PredVal(nTr,Spec,lambda2,t+1)+ouptup[4]*ADX(nA,nt+1)/100.0+ouptup[5]*MA(Ied,nMAC1,0)+ouptup[6]*MA(Ied,nMAC2,0)+ouptup[7]*MA(Ied,mMAC2,0)+ouptup[8]*MA(Stochastic(nS1,nt+1),nS1,0)/100.0+ouptup[9]*MA(Stochastic(nS2,nt+1),nS2,0)/100.0+ouptup[10]*MA(Ved,nV,0)+ouptup[11]*MA(Ved,mV,0)+ouptup[12]);
        actual.push_back(I[nt-dayahead]/I[nt]);
        i++;
    }
    return MSR(preds, actual);
};


std::vector<double> Stock::LPs(double lambda,int m, int nTr,std::vector<std::vector < double> > &Spec,int starttime,int lambda2, int nV,int mV,int nA,int nS1,int nS2,int nMAC1,int mMAC1,int nMAC2,int mMAC2,int dayahead)
{
    std::vector<double> y;
    std::vector<std::vector<double> >X;
    int t =starttime;
    int i=0;
    int j=0;
    for(j=0;j<m;j++)//cycles through all time periods
    {
        t-=20;
        std::vector<double> v;
        X.push_back(v);
        X.back().push_back(H[t]/I[t]);
        X.back().push_back(L[t]/I[t]);
        X.back().push_back(C[t]/I[t]);
        X.back().push_back(PredVal(nTr,Spec,lambda2,t));
        X.back().push_back(ADX(nA,t)/100.0);
        
        std::vector<double> Ied(mMAC2,0.0);

    for(i=0;i<mMAC2;i++)
    {
        Ied[i]=I[t+1+i]/I[t+1];
    }

    std::vector<double> Ved(mV,0.0);
    for(i=0;i<mV;i++)
    {
        Ved[i]=V[t+1+i]/V[t+1];
    }
            X.back().push_back(MA(Ied,nMAC1,0));
            X.back().push_back(MA(Ied,nMAC2,0));
            X.back().push_back(MA(Ied,mMAC2,0));

            X.back().push_back(MA(Stochastic(nS1,t),nS1,0)/100.0);
            X.back().push_back(MA(Stochastic(nS2,t),nS2,0)/100.0);
            
            X.back().push_back(MA(Ved,nV,0));
            X.back().push_back(MA(Ved,mV,0));
            
            y.push_back(I[t-dayahead]/I[t]);
    }

    const int num_dim =X[0].size();
    int num_obs = X.size();
    std::cout <<name<<" "<< X[0].size()<<" "  ;
    std::vector<double> ouptup=Lasso(num_dim,num_obs,X, y,lambda);
    for(i=0;i<X[0].size();i++)
    {
        std::cout<< ouptup[i]<< " ";
    }
    std::cout<<'\n';
    return ouptup;
};

double Stock::LPred(int t,int nTr, std::vector<std::vector < double> > &Spec,int lambda2, int nV,int mV,int nA,int nS1,int nS2,int nMAC1,int mMAC1,int nMAC2,int mMAC2,int dayahead)
{                           
    std::vector<double> Ied(mMAC2,0.0);
    int i=0;
    for(i=0;i<mMAC2;i++)
    {
        Ied[i]=I[t+1+i]/I[t+1];
    }

    std::vector<double> Ved(mV,0.0);
    for(i=0;i<mV;i++)
    {
        Ved[i]=V[t+1+i]/V[t+1];
    }
    return BetaVals[0]*H[t+1]/I[t+1]+BetaVals[1]*L[t+1]/I[t+1]+BetaVals[2]*C[t+1]/I[t+1]+BetaVals[3]*PredVal(nTr,Spec,lambda2,t+1)+BetaVals[4]*ADX(nA,t+1)/100.0+BetaVals[5]*MA(Ied,nMAC1,0)+BetaVals[6]*MA(Ied,nMAC2,0)+BetaVals[7]*MA(Ied,mMAC2,0)+BetaVals[8]*MA(Stochastic(nS1,t+1),nS1,0)/100.0+BetaVals[9]*MA(Stochastic(nS2,t+1),nS2,0)/100.0+BetaVals[10]*MA(Ved,nV,0)+BetaVals[11]*MA(Ved,mV,0);
    //+BetaVals[12]; 
};

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


std::string Stock::nameOut(void)
{
    return name;
}

Stock::Stock(std::string name1,std::vector<int> Times1, std::vector<double> I1, std::vector<double> H1, std::vector<double> L1, std::vector<double> C1, std::vector<double> V1)
{
    Times=Times1;
    name=name1;
    I=I1;
    H=H1;
    L=L1;
    C=C1;
    V=V1;
};

std::vector<double> Stock::COut(void)
{
    return C;
        
}

std::vector<double> Stock::IOut(void)
{
    return I;
        
}
std::vector<double> Stock::HOut(void)
{
    return H;
        
}
std::vector<double> Stock::LOut(void)
{
    return L;
        
}
std::vector<double> Stock::VOut(void)
{
    return V;
        
}

//need to add a function to add data to the end of the 
void Stock::PriceUpdate(double I2, double H2, double L2,double C2, double V2)
{
     std::vector<double>::iterator it;

    //it = I.begin();
    //I.insert(it, I2);
     I.push_back(I2);
    
    H.push_back(H2);
    L.push_back(L2);
    
    C.push_back(C2);
    
    V.push_back(V2);
    return;
}

//calculates moving average of vector v starting at index t for the following n elements
double MA(std::vector<double>v,int n,int t)
{
    double sum=0;
    int i=0;
    for(i=t;i<t+n;i++)
    {
        sum+=v[i];
    }
    sum=sum/n;
    return sum;
};

double EMA(std::vector<double>v,double alpha,int n,int t)
{
    double sum=v[t+n-1];
    int i=0;
    for(i=t+n-2;i>t-1;i--)
    {
        sum=alpha*v[i]+(1-alpha)*sum;
    }
    return sum;
}

std::vector<double> LassoFit(std::vector<Stock> Allstock,double lambda=1.0, int n=10,int m=25)
{
    std::vector<double> y;
    std::vector<std::vector<double> >X;
    int i =0;
    int k=0;
    int j=0;
    for(j=0;j<m;j++)//cycles through all time periods
    {
        for(k=0;k<Allstock.size();k++)//cycles through all stocks
        {
            std::vector<double> v;
            X.push_back(v);
            for(i=0;i<n;i++)//cycles through the days for each situation
            {
                X.back().push_back(Allstock[k].IOut()[j*(n+5)+i]);
                X.back().push_back(Allstock[k].HOut()[j*(n+5)+i]);
                X.back().push_back(Allstock[k].LOut()[j*(n+5)+i]);
                X.back().push_back(Allstock[k].COut()[j*(n+5)+i]);
                X.back().push_back(Allstock[k].VOut()[j*(n+5)+i]);
            }
            y.push_back(Allstock[k].IOut()[j*(n+5)]);
        }
    }
    const int num_dim =X[0].size();
    int num_obs = X.size();
    return Lasso(num_dim,num_obs,X, y,lambda);
    
}

/*std::vector<double> Stock::Trend(int n,int lambda,int t)
{
    
    int i=0;
        int j=0;

    
    //for hodrick prescott filter estimate
        
    std::vector<std::vector<double> > Spec(n, std::vector<double>(n,0.0));  
    std::vector<std::vector<double> > K(n, std::vector<double>(n,0.0)); 
    //std::vector<std::vector<double> > Spec(HistLength()-t, std::vector<double>(HistLength()-t,0.0));
    
    //std::vector<std::vector<double> > K(HistLength()-t,std::vector<double>(HistLength()-t,0.0));
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
               
    std::vector<double> tau(n,0.0);
    i=0;
    j=t;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            tau[i]+=Spec[i][j]*log(C[t+j]);
        }
        tau[i]=exp(tau[i]);
        //std::cout << tau[i]<< " " ;
    }
    //std::cout <<'\n';
    return tau;
      
};*/


std::vector<double> Stock::Trend(int n,std::vector<std::vector<double> >Spec, int lambda,int t)
{
    
    int i=0;
    int j=t;

    
    //for hodrick prescott filter estimate
    //The matrix Spec is the desired matrix to use, and has been calculated separately for speed
               
    std::vector<double> tau(n,0.0);//vector with trend component
    
    //create normalised values of I to use. Normalise based on value of I[t]
    //do i need to scale by standard deviation as well?
    std::vector<double> Ied(n,0.0);
    for(i=0;i<n;i++)
    {
        Ied[i]=I[t+i]/I[t];
    }
    for(i=0;i<n;i++)
    {
        //multiplies matrix by the vector of the previous n days closing data to obtain a smoothed trend line through the stock prices
        for(j=0;j<n;j++)
        {
            tau[i]+=Spec[i][j]*Ied[j];
        }

    }

    return tau;
      
};

double Stock::PredVal(int n,std::vector<std::vector<double> >Spec, int lambda,int t)
{
    std::vector<double> tau=Trend(n,Spec, lambda, t);
    double sum=0;
    int i=0;
    for(i=0;i<n-1;i++)
    {
        sum+=tau[i]-tau[i+1];
    }
    sum/=((double)n-(double)1);
    return (n*sum)/tau[0];
}

int Stock::TrendTest(int n,std::vector<std::vector<double> >Spec,int lambda,int t)
{

    int i=0;

    if(Stock::Trend(n,Spec,lambda,t)[0]-Stock::Trend(n,Spec,lambda,t)[1]>0)
    {
        
        for(i=1;i<n-1;i++)
        {

            if(Stock::Trend(n,Spec,lambda,t)[i]-Stock::Trend(n,Spec,lambda,t)[i+1]<=0)
            {

                return 0;
            }  
            
        }

        return 1;
    }
    if(Stock::Trend(n,Spec,lambda,t)[0]-Stock::Trend(n,Spec,lambda,t)[1]<0)
    {
        for(i=1;i<n-1;i++)
        {

            if(Stock::Trend(n,Spec,lambda,t)[i]-Stock::Trend(n,Spec,lambda,t)[i+1]>=0)
            {

                return 0;
            }
            
        }

        return -1;
    }
    if(Stock::Trend(n,Spec,lambda,t)[0]-Stock::Trend(n,Spec,lambda,t)[1]==0)
    {
        return 0;
    }
}


//produces volume test for stock, with parameters n and m
//defaults should be n=2 and m =14
bool Stock::VolTest(int n,int m,int t)
{
    bool outcome=false;
    //if(EMA(V,2/(n+1),n,t)>EMA(V,2/(m+1),m,t))
    if(MA(V,n,t)>MA(V,m,t))
    {
        outcome=true;
    }
    return outcome;
};

std::vector<double> Stock::DMplus(int n,int t)
{
    std::vector<double> out(n,0.0);
    int i=0;
    for(i=0;i<n;i++)
    {
        if(H[t+i]-H[t+i+1]>=L[t+i+1]-L[t+i])
        {
            out[i]=std::max(H[t+i]-H[t+i+1],(double) 0);
        }
    }
    return out;
};
std::vector<double> Stock::DMminus(int n,int t)
{
    std::vector<double> out(n,0.0);
    int i=0;
    for(i=0;i<n;i++)
    {
        if(H[t+i]-H[t+i+1]<=L[t+1+i]-L[t+i])
        {
            out[i]=std::max(L[t+i]-L[t+i+1],(double) 0);
        }

    }
    return out;
};
std::vector<double> Stock::TR(int n,int t)
{
    std::vector<double> out(n,0.0);
    int i=0;
    for(i=0;i<n;i++)
    {
        double a =H[t+i]-L[t+i];
        double b=fabs(H[t+i]-C[t+i+1]);
        double c =fabs(L[t+i]-C[t+i+1]);
        out[i]=std::max(std::max(a,b),c);
    }
    return out;
};
std::vector<double> Stock::DIplus(int n,int t)
{
    std::vector<double> out(n,0.0);
    int i=0;
    for(i=0;i<n;i++)
    {
        double a =Stock::DMplus(n,t)[i];
        double b=Stock::TR(n,t)[i];
        out[i]=a/b;
    }
    
    //std::cout <<out.size()<<'\n';
    return out;
};
std::vector<double> Stock::DIminus(int n,int t)
{
    std::vector<double> out(n,0.0);
    int i=0;
    for(i=0;i<n;i++)
    {
        out[i]=Stock::DMminus(n,t)[i]/Stock::TR(n,t)[i];

    }
    return out;
};
//usual parameter is n=14
//true means strong trend,
//false means weak trend
double Stock::ADX(int n,int t)
{
    std::vector<double> val(n,0.0);
    int i=0;
    for(i=0;i<n;i++)
    {

        //val[i]=(EMA(Stock::DIplus(2*n,t),2/(n+1),n,i)-EMA(Stock::DIminus(2*n,t),2/(n+1),n,i))/(EMA(Stock::DIplus(2*n,t),2/(n+1),n,i)+EMA(Stock::DIminus(2*n,t),2/(n+1),n,i));
        val[i]=(MA(Stock::DIplus(2*n,t),n,i)-MA(Stock::DIminus(2*n,t),n,i))/(MA(Stock::DIplus(2*n,t),n,i)+MA(Stock::DIminus(2*n,t),n,i));
    }
    //return 100*EMA(val,2/(n+1),n,0);
    return 100*MA(val,n,0);
}
bool Stock::ADXTest(int n,int t)
{
    std::vector<double> val(n,0.0);
    int i=0;
    for(i=0;i<n;i++)
    {
            //val[i]=(EMA(Stock::DIplus(2*n,t),2/(n+1),n,i)-EMA(Stock::DIminus(2*n,t),2/(n+1),n,i))/(EMA(Stock::DIplus(2*n,t),2/(n+1),n,i)+EMA(Stock::DIminus(2*n,t),2/(n+1),n,i));
        val[i]=(MA(Stock::DIplus(2*n,t),n,i)-MA(Stock::DIminus(2*n,t),n,i))/(MA(Stock::DIplus(2*n,t),n,i)+MA(Stock::DIminus(2*n,t),n,i));
    }
    double ADX=100*MA(val,n,0);
    //double ADX=100*EMA(val,2/(n+1),n,0);
    if(ADX>25)
    {
        return true;
    }
    if(ADX<=25)
    {
        return false;
    }
    return false;
};

std::vector<double> Stock::Stochastic(int n,int t)
{
    std::vector<double> K;
    int i=0;
    for(i=t;i<t+5*n;i++)
    {
        double Ln=L[i];
        double Hn=H[i];
        int j=0;
        for( j=0;j<n;j++)
        {
            Ln=std::min(Ln,L[i+j]);
            Hn=std::max(Hn,H[i+j]);
        }
        K.push_back(100*(C[i]-Ln)/(Hn-Ln));
    }
    return K;
};

int Stock::StochTest(int n,int t)
{
    std::vector<double> FD;
    int i=0;
    
    for(i=0;i<6;i++)
    {
        FD.push_back(MA(Stock::Stochastic(n,t),3,i));
        //FD.push_back(EMA(Stock::Stochastic(n,t),0.5,3,i));
    }
    
    double D=MA(FD,3,0);
    double D2=MA(FD,3,1);
    //double D=EMA(FD,0.5,3,0);
    //double D2=EMA(FD,0.5,3,1);
    
    if((Stock::Stochastic(n,t)[1]<D2)&&(D2<20))
    {
        if((D<Stock::Stochastic(n,t)[0])&&(Stock::Stochastic(n,t)[0]<20))
        {
            return 1;
        }
    }
    
    if((Stock::Stochastic(n,t)[1]>D2)&&(D2>80))
    {
        if((D>Stock::Stochastic(n,t)[0])&&(Stock::Stochastic(n,t)[0]>80))
        {
            return -1;
        }   
    }
    else
    {
        return 0;
    }
};


//outcome, 1 is buy, 0 is hold, -1 is sell
int Stock::MACTest(int n, int m,int t)
{
    int total =0;
    //check if there is a crossover for each price stream
    //opening price
    if(MA(I,n,t+1)<MA(I,m,t+1))
    {
        if(MA(I,n,t)>MA(I,m,t))
        {
            total+= 1;
        }
    }
    if(MA(I,n,t+1)>MA(I,m,t+1))
    {
        if(MA(I,n,t)<MA(I,m,t))
        {
            total =total -1;
        }
    }
    
    /*//high price
    if(MA(H,n,t+1)<MA(H,m,t+1))
    {
        if(MA(H,n,t)>MA(H,m,t))
        {
            total+= 1;
        }
    }
    if(MA(H,n,t+1)>MA(H,m,t+1))
    {
        if(MA(H,n,t)<MA(H,m,t))
        {
            total =total -1;
        }
    }
    //low price
        if(MA(L,n,t+1)<MA(L,m,t+1))
    {
        if(MA(L,n,t)>MA(L,m,t))
        {
            total+= 1;
        }
    }
    if(MA(L,n,t+1)>MA(L,m,t+1))
    {
        if(MA(L,n,t)<MA(L,m,t))
        {
            total =total -1;
        }
    }*/
    
    //close price
    if(MA(C,n,t+1)<MA(C,m,t+1))
    {
        if(MA(C,n,t)>MA(C,m,t))
        {
            total+= 1;
        }
    }
    if(MA(C,n,t+1)>MA(C,m,t+1))
    {
        if(MA(C,n,t)<MA(C,m,t))
        {
            total =total -1;
        }
    }
    
    return total;
    
    //using EMA instead of SMA
    /*
    if(EMA(I,2/(n+1)n,t+1)<EMA(I,2/(m+1),m,t+1))
    {
        if(EMA(I,2/(n+1),n,t)>EMA(I,2/(m+1),m,t))
        {
            total+= 1;
        }
    }
    if(EMA(I,2/(n+1),n,t+1)>EMA(I,2/(m+1),m,t+1))
    {
        if(EMA(I,2/(n+1),n,t)<EMA(I,2/(m+1),m,t))
        {
            total =total -1;
        }
    }
     * */
    
};


int Stock::CandTest(int n,std::vector<std::vector<double> >Spec,int lambda,int t)
{
    
    int i=0;
    //k is the maximum number of previous days one needs
    //currently only require 3
    int k=3;
    //First Create Various parameters to describe the candlestick patterns
    
    std::vector<bool> Colour(k,0);
    for(i=0;i<k;i++)
    {
        if(C[t+i]>I[t+i])
        {
            Colour[i]=1;
        }
    }
    
    std::vector<double> Range(k,0.0);
    for(i=0;i<k;i++)
    {
        Range[i]=H[t+i]-L[t+i];
    }
    
    std::vector<double> Body(k,0.0);
    for(i=0;i<k;i++)
    {
        double a=fabs(C[t+i]-I[t+i]);
        Body[i]=a/Range[i];
    }
    std::vector<double> UShad(k,0.0);
    for(i=0;i<k;i++)
    {
        if(Colour[i]==1)
        {
            UShad[i]=(H[t+i]-C[t+i])/Range[i];
        }
        if(Colour[i]==0)
        {
            UShad[i]=(H[t+i]-I[t+i])/Range[i];
        }
    }
    
    std::vector<double> LShad(k,0.0);
    for(i=0;i<k;i++)
    {
        if(Colour[i]==1)
        {
            UShad[i]=(I[t+i]-L[t+i])/Range[i];
        }
        if(Colour[i]==0)
        {
            UShad[i]=(C[t+i]-L[t+i])/Range[i];
        }
    }
    
    //we now go through various patterns to check if one has achieved it.
    
    
    //First Pattern is the Hammer
    int a =TrendTest(n,Spec,lambda,t);
    if(a==-1)
    {
        if((Colour[0]==1)&&(LShad[0]>=2*Body[0])&&(UShad[0]<=0.05)&&(Body[0]>=0.1))
        {
            return 1;
        }
    }
    //Second Pattern is Bullish Engulfing pattern
    if(a==-1)
    {
        if((Colour[1]==0)&&(Colour[0]==1)&&(I[t]<=C[t+1])&&(C[t]>=I[t-1]))
        {
            return 1;
        }
    }
    
    //Hanging Man
    if(a==1)
    {
        if((Colour[0]==0)&&(LShad[0]>=2*Body[0])&&(UShad[0]<=0.05)&&(Body[0]>=0.1))
        {
            return -1;
        }
    }
    //Bearish Engulfing
    if(a==-1)
    {
        if((Colour[1]==1)&&(Colour[0]==0)&&(I[t]>=C[t+1])&&(C[t]<=I[t-1]))
        {
            return -1;
        }
    }
    //The default case is that no patterns have been seen and so 
    return 0;
};

//outputs integer
//1 means Buy
//-1 means sell
// 0 means hold
int Stock::BuySellDec(int t,std::vector<std::vector<double> >Spec,int nTr,int lambda2, int nV,int mV,int nA, int nS1, int nS2, int nMAC1,int mMAC1,int nMAC2,int mMAC2,int dayahead)
{
    int result=0;
    double value = LPred(t, nTr, Spec, lambda2, nV, mV, nA, nS1, nS2, nMAC1, mMAC1, nMAC2, mMAC2,dayahead);
    //std::cout << name << " " << value<<'\n';
    if(value >1.05)
    {
        result =1;
    }
    if(value<1)
    {
        result =-1;
    }
    if((value >=1)&&(value <=1.05))
    {
        result =0;
    }
    return result;
}

/*
int Stock::BuySellDec(int t,std::vector<std::vector<double> >Spec,int nTr,int lambda, int nV,int mV,int nA, int nS1, int nS2, int nMAC1,int mMAC1,int nMAC2,int mMAC2,int u1,int d1,int u2,int d2,int u3,int d3,double Tu,double frac)
{
    
    //at present time t, the decision to buy or sell is based on the previous n days stock values, when n ranges from the day before t ...
    //therefore for this algorithm require time to range from t+1...
    int tt= t+1;
    int a =Stock::TrendTest(nTr,Spec,lambda,tt);
    bool b=Stock::VolTest(nV,mV,tt);
    bool c=Stock::ADXTest(nA,tt);
    int d=Stock::StochTest(nS1,tt)+Stock::StochTest(nS2,tt) + Stock::MACTest(nMAC1, mMAC1,tt) + Stock::MACTest(nMAC2, mMAC2,tt)+Stock::CandTest(nTr,Spec,lambda,tt);
    double f=PredVal(nTr,Spec,lambda,tt);
    
    //test for Uptrend in last 20 days
    if(a==1)
    {
        //if uptrend, then decide whether time is right to buy or not
        
        //Tests of b and c are strength of trend
        //first case is strong strend
        if((b==true)&&(c==true))
        {
            //so long as enough terms in d predicts now is a good time to buy, then buy
            if(d>=u1)
            {
                //if prediction of return is sufficiently high, then buy
                if(f>=Tu)
                {
                    return 1;
                }
            }
            
            //if predict now is a good time to sell, then sell
            if(d<=d1)
            {
                    
                return -1;
                    
            }
                
        }
        //semi strong trend when one says strong and other doesnt
        if(((b==true)&&(c==false))||((b==false)&&(c==true)))
        {
                if(d>=u2)
                {
                    if(f>=Tu)
                    {
                        return 1;
                    }
                }
                if(d<=d2)
                {

                    return -1;
                    
                }
        }
            
        
        if((b==false)&&(c==false))
        {
            if(d>=u3)
            {
                if(f>=Tu)
                {
                    return 1;
                }
            }
            
            if(d<=d3)
            {
                return -1;
                    
            }
        }    
    }
    
    //if trend is not up, then sell the stock
    if(a==0)
    {
        return -1;

    }
    if(a==-1)
    {
        return -1;
    } 
    
    //if none of the above Buy/Sell decision criteria have been met, then we do nothing.
    //in case the options above do not encompass all choices
    return 0;
};*/

int Stock::HistLength(void)
{
    return Times.size();
}

double StockLatHigh(std::string name, std::vector<Stock> &St, int time)
{
    double val=0;
    int k =0;
    for(k=0;k<St.size();k++)
    {
        
        if(St[k].nameOut().find(name) != std::string::npos)
        {
            
            val=St[k].HOut()[time];
        }
        
    }
    
    return val;
    
}
double StockLatLow(std::string name, std::vector<Stock> &St, int time)
{
    double val=0;
    int k =0;
    for(k=0;k<St.size();k++)
    {
        
        if(St[k].nameOut().find(name) != std::string::npos)
        {
            
            val=St[k].LOut()[time];
        }
        
    }
    
    return val;
    
}
double StockLatOp(std::string name, std::vector<Stock> &St,int time)
{
    double val=0;
    int k =0;
    for(k=0;k<St.size();k++)
    {
        
        if(St[k].nameOut().find(name) != std::string::npos)
        {
            
            val=St[k].IOut()[time];
        }
        
    }
    
    return val;
    
}
double StockLatCl(std::string name, std::vector<Stock> &St,int time)
{
    double val=0;
    int k =0;
    for(k=0;k<St.size();k++)
    {
        
        if(St[k].nameOut().find(name) != std::string::npos)
        {
            
            val=St[k].COut()[time];
        }
        
    }
    
    return val;
    
}
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//for outputting and inputting data

//still need a routine for inputting data downloaded from IB API into the class I have created
 int NStockInput(std::vector<Stock> &Input)
{
     int i=0;
    
    std::ifstream seclist;
    seclist.open("Data//seclist.txt",std::ios::in);
    if (!seclist) 
    {//ensures that the file can indeed be opened.
            std::cout << "==============================================================="<<'\n';
            std::cout << " Unable to open file seclist.txt" << '\n';
            std::cout << " Therefore we are closing the program." << '\n';
            std::cout << " " <<'\n';
            std::cout << " The user should be aware that there must be a file called "<<'\n';
            std::cout << " seclist.txt in the same directory as the program for " <<'\n';
            std::cout << " the program to work. Maybe create a file and try again?"<<'\n';
            std::cout << "==============================================================="<<'\n';
        return 0;
    }
    
    std::vector<std::string> securitylist;
    while(seclist.good())
    {
        std::string line; //creates empty string
        getline(seclist, line); //reads file into the empty string
        if(line.length()>=2)
        {
            securitylist.push_back(line);   
        }
    }
    seclist.close();
    std::ifstream input;
    for(i=0;i<securitylist.size();i++)
    {
        std::string filename="Data//";
        filename.append(securitylist[i]);
        filename.append(".L.csv");
            
        input.open(filename.c_str(),std::ios::in);
        if (!input) 
        {//ensures that the file can indeed be opened.
            
            std::cout << "==============================================================="<<'\n';
            std::cout <<  "Cant open "<< filename.c_str() << '\n';
            std::cout << "==============================================================="<<'\n';
        }
        std::vector<int> Tim;
        std::vector<double> Start;
        std::vector<double> High;
        std::vector<double> Low;
        std::vector<double> Close;
        std::vector<double> Volume;
        
        while(input.good())
        {
            std::string line;
            getline(input,line);
            
            if(line.size()==0)
            {break;}
            std::size_t found = line.find(",");
            std::string dateval= line.substr(0,found);
            line=line.erase(0,found+1);
            found =dateval.find("/");
            std::string temp=dateval.substr(0,found);
            dateval=dateval.erase(0,found+1);
            int month = stoi(temp);
            
            found = dateval.find("/");
            temp=dateval.substr(0,found);
            dateval=dateval.erase(0,found+1);
            int day =stoi(temp);
            
            temp=dateval;
            int year = stoi(temp);
            
            Day Date(year, month,day);
            //std::cout << year<< " " << month << " " << day<< " " << Date.count<<'\n';
            std::vector<int>::iterator it2;
            it2=Tim.begin();
            it2=Tim.insert(it2, Date-Reference);
            
            std::vector<double> vect;
            while(line.size()>0)
            {
                found =line.find(",");
                temp = line.substr(0,found);
                double num = stof(temp);
                vect.push_back(num);
                if(found+1<=line.size())
                {
                    line.erase(0,found+1);
                }
                    
                if(vect.size()==5)
                {
                    break;
                }           
            }
            if(vect.size()==5)
            {
                if(vect[4]!=0)
                {
                     std::vector<double>::iterator it;
                    it = Start.begin();
                    it=Start.insert(it,vect[0]);
                    
                    it =High.begin();
                    it=High.insert(it,vect[1]);
                    
                    it=Low.begin();
                    it=Low.insert(it,vect[2]);
                    
                    it=Close.begin();
                    it=Close.insert(it,vect[3]);
                    
                    it=Volume.begin();
                    it=Volume.insert(it,vect[4]);
                }
                if(vect[4]==0)
                {
                    Tim.erase(Tim.begin(),Tim.begin()+1);
                }
            }
            
        }
        quickSort4(Tim, Start, High, Low, Close, Volume, 0, Tim.size());
        std::string B=securitylist[i];
        Stock New(B,Tim, Start, High, Low, Close, Volume);
        Input.push_back(New);
        input.close();
        
    }
    return Input.size();
}
 int StockInput(std::vector<Stock> &Input)
{
    
    int i=0;
    
    std::ifstream seclist;
    seclist.open("Data//seclist.txt",std::ios::in);
    if (!seclist) 
    {//ensures that the file can indeed be opened.
            std::cout << "==============================================================="<<'\n';
            std::cout << " Unable to open file seclist.txt" << '\n';
            std::cout << " Therefore we are closing the program." << '\n';
            std::cout << " " <<'\n';
            std::cout << " The user should be aware that there must be a file called "<<'\n';
            std::cout << " seclist.txt in the same directory as the program for " <<'\n';
            std::cout << " the program to work. Maybe create a file and try again?"<<'\n';
            std::cout << "==============================================================="<<'\n';
        return 0;
    }
    
    std::vector<std::string> securitylist;
    while(seclist.good())
    {
        std::string line; //creates empty string
        getline(seclist, line); //reads file into the empty string
        if(line.length()>=2)
        {
            securitylist.push_back(line);   
        }
    }
    seclist.close();
        
    std::ifstream input;
    for(i=0;i<securitylist.size();i++)
    {
        std::string filename="Data//";
        filename.append(securitylist[i]);
        filename.append(".txt");
            
        input.open(filename.c_str(),std::ios::in);
        if (!input) 
        {//ensures that the file can indeed be opened.
            
            std::cout << "==============================================================="<<'\n';
            std::cout <<  "Cant open "<< filename.c_str() << '\n';
            std::cout << "Creating Such a file"<<'\n';
            std::cout << "==============================================================="<<'\n';
            std::ofstream newfile;
            newfile.open(filename.c_str(),std::ios::out);
            newfile.close();

        }
        std::vector<int> Tim;
        std::vector<double> Start;
        std::vector<double> High;
        std::vector<double> Low;
        std::vector<double> Close;
        std::vector<double> Volume;
        
        while(input.good())
        {
            std::string line;
            getline(input,line);
            std::stringstream liner(line);
            double num=0;
            std::vector<double> info;
            while(liner >>num)
            {
                
                info.push_back(num);
                
            }
                
            if(info.size()==8)
            {
                Day Date(info[0],info[1],info[2]);
                Tim.push_back(Date-Reference);
                Start.push_back(info[3]);

                High.push_back(info[4]);
                Low.push_back(info[5]);
                Close.push_back(info[6]);
                Volume.push_back(info[7]);
            }
            
            if(info.size()==5)
            {
                Start.push_back(info[0]);

                High.push_back(info[1]);
                Low.push_back(info[2]);
                Close.push_back(info[3]);
                Volume.push_back(info[4]);
            }
                    
        }
        
        std::string B=securitylist[i];
        Stock New(B,Tim, Start, High, Low, Close, Volume);
        Input.push_back(New);
        input.close();
    }

     
    return Input.size();
};

int WriteStockList(std::vector<Stock> Input)
{
    int i=0;
    std::ofstream output;
    for(i=0;i<Input.size();i++)
    {
        std::string filename="Data//";
        filename.append(Input[i].nameOut());
        filename.append(".txt");
            
        output.open(filename.c_str(),std::ios::out);
        
        int j=0;
        for(j=0;j<Input[i].COut().size();j++)
        {
            output <<Dayinv (Input[i].Times[j])[0]<<" " <<Dayinv (Input[i].Times[j])[1] <<" " << Dayinv (Input[i].Times[j])[2]<<" "  <<Input[i].IOut()[j]<< " " << Input[i].HOut()[j]<< " " << Input[i].LOut()[j]<< " " << Input[i].COut()[j]<< " " << Input[i].VOut()[j]<<'\n';
        }
        output.close();
    }
    return 0;
}
/////////////////////////////////////////////////////////////////////////////////
//produces a list of all possible stocks with their buy sell data
void StockOutput(std::vector<Stock> &Universe,int nTr,std::vector<std::vector<double> > &Spec,int lambda, int nV,int mV,int nA,int nS1,int nS2,int nMAC1,int mMAC1,int nMAC2,int mMAC2, int t,double FC,double Cost,double frac1,int dayahead)
{
    int i=0;

    std::ofstream Sellfile;
    std::string filename="Outputs//";
    filename.append("Sell");
    filename+=std::to_string(Cost);
    filename.append("-");
    filename+=std::to_string(frac1);
    filename.append("-");
    filename+=std::to_string(dayahead);
    filename.append(".txt");

    Sellfile.open (filename.c_str(), std::ios::out);
    
    
    std::vector<std::string> Buylist;
    std::vector<double> Strength;
    while(i<Universe.size())
    {
        int a =Universe[i].BuySellDec(t,Spec, nTr, lambda,  nV,mV, nA, nS1,nS2, nMAC1, mMAC1, nMAC2, mMAC2,dayahead);
        if(a==-1)
        {
            
            Sellfile<< Universe[i].nameOut()<<'\n';
            //std::cout <<"SELL " << Universe[i].nameOut()<<'\n';

        } else{
            if(a==1)
            {
                Buylist.push_back(Universe[i].nameOut());
                //for Strength have multiple choices
                Strength.push_back(Universe[i].ADX(nA, t));
                //Strength.push_back(Universe[i].ADXTest(nA,t)+Universe[i].VolTest(nV,mV,t));
                //Strength.push_back(Universe[i].PredVal(nTr,Spec,lambda,t));

            }
        }
        i++;
    }
    
    Sellfile.close();
    
    quickSort(Strength,Buylist,0,Strength.size());
    
    //write buy list to file
    std::ofstream Buyfile;
    filename="Outputs//";
    filename.append("Buy");
    filename+=std::to_string(Cost);
    filename.append("-");
    filename+=std::to_string(frac1);
    filename.append("-");
    filename+=std::to_string(dayahead);
    filename.append(".txt");
    
    Buyfile.open (filename.c_str(), std::ios::out);
    if (!Buyfile) 
    {//ensures that the file can indeed be opened.
        std::cout << "==============================================================="<<'\n';
        std::cout << " Unable to open file Buy.txt" << '\n';
        std::cout << "==============================================================="<<'\n';
    }
    for(i=0;i<Buylist.size();i++)
    {
        Buyfile<< Buylist[i]<<'\n';
         //std::cout <<"BUY  " << Buylist[i]<<'\n';
    }
        

    Buyfile.close();

}

//produces a file with a collection of stocks to buy
void StockOutputBuy(std::vector<Stock> &Universe,int nTr,std::vector<std::vector<double> > &Spec,int lambda, int nV,int mV,int nA,int nS1,int nS2,int nMAC1,int mMAC1,int nMAC2,int mMAC2,int t,double FC,double Cost,double frac1,int dayahead)
{
    int i=0;
    
    std::vector<std::string> Buylist;
    std::vector<double> Strength;
    //std::cout << Universe.size()<<'\n';
    for(i=0;i<Universe.size();i++)
    {
        /*//Case 1 is order of whatever Universe is ordered buy
        if(Universe[i].BuySellDec(t,Spec, nTr, lambda,  nV,mV, nA, nS, nMAC1, mMAC1, nMAC2, mMAC2,u11,d11, u12, d12, u13, d13, Tu1)==1)
        {
            Outputfile<< Universe[i].nameOut()<<'\n';
        }*/
        
        
        //Case 2 is order based on strength of something. Free to choose values
        if(Universe[i].BuySellDec(t,Spec, nTr, lambda,  nV,mV, nA, nS1,nS2, nMAC1, mMAC1, nMAC2, mMAC2, dayahead)==1)
        {
            Buylist.push_back(Universe[i].nameOut());
            //for Strength have multiple choices
            Strength.push_back(Universe[i].ADX(nA, t));
            //Strength.push_back(Universe[i].ADXTest(nA,t)+Universe[i].VolTest(nV,mV,t));
            //Strength.push_back(Universe[i].PredVal(nTr,Spec,lambda,t));

        }
    }
    
    //one then sorts based on the values in Strength
    quickSort(Strength,Buylist,0,Strength.size());
    
    std::ofstream Outputfile;
    std::string filename="Outputs//";
    filename.append("Buy");
    filename+=std::to_string(Cost);
    filename.append("-");
    filename+=std::to_string(frac1);
    filename.append(".txt");
    Outputfile.open (filename.c_str(), std::ios::out);
    if (!Outputfile) 
    {//ensures that the file can indeed be opened.
        std::cout << "==============================================================="<<'\n';
        std::cout << " Unable to open file Buy.txt" << '\n';
        std::cout << "==============================================================="<<'\n';
    }
    
    for(i=0;i<Buylist.size();i++)
    {
        Outputfile<< Buylist[i]<<'\n';
    }
    Outputfile.close();
    return;
}


//produces a file with a list of stocks to sell
void StockOutputSell(std::vector<Stock> &Universe,int nTr,std::vector<std::vector<double> > &Spec,int lambda, int nV,int mV,int nA,int nS1,int nS2,int nMAC1,int mMAC1,int nMAC2,int mMAC2,int t,double FC,double Cost,double frac1,int dayahead)
{
    int i=0;
    
    std::ofstream Outputfile;
    std::string filename="Outputs//";
    filename.append("Sell");
    filename+=std::to_string(Cost);
    filename.append("-");
    filename+=std::to_string(frac1);
    filename.append(".txt");
    Outputfile.open (filename.c_str(), std::ios::out);
    for(i=0;i<Universe.size();i++)
    {
        if(Universe[i].BuySellDec(t,Spec, nTr, lambda, nV,mV, nA, nS1,nS2,nMAC1, mMAC1, nMAC2, mMAC2, dayahead)==-1)
        {
            
            Outputfile<< Universe[i].nameOut()<<'\n';

        }
        
    }
    Outputfile.close();
    return;
}
