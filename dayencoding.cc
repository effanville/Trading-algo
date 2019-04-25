#include <vector>
#include <iostream>
#include "dayencoding.h"


//creator for a day
int Day::month_days[] = {0,31,59,90,120,151,181,212,243,273,304,334};

//Create a reference date for which all other dates are compared to
Day Reference(2000,1,1);

//default constructor for an element of Day
Day::Day (int day_year,
        int day_month,
        int day_days){
        // calculate number of leap years.
        int leapyears    = day_year / 4;
        if (day_year % 4 == 0 && day_month < 3)
        {
            leapyears   --;
        }
        // convert year/month/day into a day count
        count    = day_year * 365 + month_days[day_month-1] + day_days + leapyears;
}

//provide inverse for the day structure, given that the input is the number of days from 1 1 2000
//works by calculating the number of years that have passed since 2000,
//and then by calculating the number of months in that year, based upon 
//knowing number of days in each month
std::vector<int> Dayinv (int d)
{
    std::vector<int> opt;
    //provide estimate for the number of years, up to calculating leap years
    int ny=0;
    int ly=0;
    //iteratively add to number of years until reach desired one that has max number of days below target d
    while(365*ny+ly<=d)
    {
        if(ny%4==0)
        {
            ly++;
        }
        ny++;
    }
    //the above while loop adds 1 onto ly and ny beyond what is needed, so must remove this
    ny--;
    if(ny%4==0)
    {
        ly--;
    }
    opt.push_back(2000+ny);//input the year as the first element of the output vector
    int rem = d-365*ny-ly;//find number of days seen in current year
    int nm=0;
    
    int i=0;
    int rem2;
    int year;
    year =2000+ny;
    if (year % 4 == 0)//number of days for end of month is different in a leap year
    {
        int md[]={0,31,60,91,121,152,182,213,244,274,305,335};
        for(i=0; i<12; i++)//cycle through months until find the current one
        {
            if( rem -md[i]>=0)
            {
                nm=i+1;
                rem2=1+rem-md[i];//find the day of the month
            }
        }
    }
    
    if(year % 4!=0)//if not in a leap year
    {
        int md[]={0,31,59,90,120,151,181,212,243,273,304,334};
        for(i=0; i<12; i++)
        {
            if( rem -md[i]>=0)
            {
                nm=i+1;
                rem2=1+rem-md[i];
            }
        }
    }
    opt.push_back(nm);//inputs number of months
    opt.push_back(rem2);//inputs number of months
    return opt;
}
