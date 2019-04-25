#include <vector>
#include <iostream>
#include <string>

////////////////////////////////////////////////
//to perform quicksort algorithm on data, sorting via time A and ensuring value B is updated accordingly

//splitting of vectors into parts larger and smaller than the selected pivot value (taken to be the value A[p]
//sorts in order largest to smallest
int partition4(std::vector<int>& A,std::vector<double>& B,std::vector<double>& C,std::vector<double>& D,std::vector<double>& E,std::vector<double>& F, int p,int q)
{
    
    
    int x= A[p];
    int i=p;
    int j;

    for(j=p+1; j<q; j++)
    {
        if(A[j]>=x)
        {
            i=i+1;
            std::swap(A[i],A[j]);
            std::swap(B[i],B[j]);
            std::swap(C[i],C[j]);
            std::swap(D[i],D[j]);
            std::swap(E[i],E[j]);
            std::swap(F[i],F[j]);
        }

    }

    std::swap(A[i],A[p]);
        std::swap(B[i],B[p]);
        std::swap(C[i],C[p]);
        std::swap(D[i],D[p]);
        std::swap(E[i],E[p]);
                std::swap(F[i],F[p]);
    return i;
}

//main iterative algorithm for quicksort. involves first partitioning the vectors, and then quicksorting the vectors above and below the pivot
//default pivot is the first element of A
void quickSort4(std::vector<int>& A,std::vector<double>& B,std::vector<double>& C,std::vector<double>& D,std::vector<double>& E,std::vector<double>& F, int p,int q)
{
    if(A.size()!=B.size())
    {
        std::cout << " Incompatible vectors"<<'\n';
        return;
    }
    int r;
    if(p<q)
    {
        r=partition4(A,B,C,D,E,F, p,q);
        quickSort4(A,B,C,D,E,F,p,r);  
        quickSort4(A,B,C,D,E,F,r+1,q);
    }
}

//splitting of vectors into parts larger and smaller than the selected pivot value (taken to be the value A[p]
int partition(std::vector<double>& A,std::vector<std::string>& B, int p,int q)
{
    
    
    int x= A[p];
    int i=p;
    int j;

    for(j=p+1; j<q; j++)
    {
        if(A[j]>=x)
        {
            i=i+1;
            std::swap(A[i],A[j]);
            std::swap(B[i],B[j]);
        }

    }

    std::swap(A[i],A[p]);
        std::swap(B[i],B[p]);
    return i;
}

//main iterative algorithm for quicksort. involves first partitioning the vectors, and then quicksorting the vectors above and below the pivot
//default pivot is the first element of A
void quickSort(std::vector<double>& A,std::vector<std::string>& B, int p,int q)
{
    if(A.size()!=B.size())
    {
        std::cout << " Incompatible vectors"<<'\n';
        return;
    }
    int r;
    if(p<q)
    {
        r=partition(A,B, p,q);
        quickSort(A,B,p,r);  
        quickSort(A,B,r+1,q);
    }
}
