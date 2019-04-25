#pragma once

#ifndef quicksort_H
#define quicksort_H

int partition(std::vector<double>& A,std::vector<std::string>& B, int p,int q);

void quickSort(std::vector<double>& A,std::vector<std::string>& B, int p,int q);

int partition4(std::vector<int>& A,std::vector<double>& B,std::vector<double>& C,std::vector<double>& D,std::vector<double>& E,std::vector<double>& F, int p,int q);

void quickSort4(std::vector<int>& A,std::vector<double>& B,std::vector<double>& C,std::vector<double>& D,std::vector<double>& E,std::vector<double>& F, int p,int q);

#endif
