#pragma once

#ifndef dayencoding_H
#define dayencoding_H


//////////////////////////////////////////////////////////
//structure for storing a date as an integer. 

struct Day
{
    int count; //only variable required is number of days from some arbitrary start
    //input operator that takes an input stream and a member of Day
    friend std::istream& operator>>(std::istream& s, Day& d)
    {
        int day_year;
        int day_month;
        int day_days;
        // calculate number of leap years.
        int leapyears    = day_year / 4;
        if (day_year % 4 == 0 && day_month < 3)
        {
            // If this is a leap year and we have not passed Feburary then it does not count
            leapyears   --;
        }
        // convert year/month/day into a day countk
        d.count    = day_year * 365 + month_days[day_month-1] + day_days + leapyears;

        // return stream for chaining
        return s;
    }
    friend int operator-(Day const& lhs, Day const& rhs)
    {
        // subtraction gives you the difference between two Days objects.
        return lhs.count - rhs.count;
    }
    Day( int day_year, int day_month, int day_days);

    static int month_days[];
};

//provide inverse for the day structure, given that the input is the number of days from 1 1 2000
//works by calculating the number of years that have passed since 2000,
//and then by calculating the number of months in that year, based upon 
//knowing number of days in each month
std::vector<int> Dayinv (int d);

#endif
