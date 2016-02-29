//
//  scrub.h
//  9898BDiFHW1
//
//  Created by Xiaobo on 2/20/16.
//  Copyright © 2016 Xiaobo. All rights reserved.
//

#ifndef scrub_HPP
#define scrub_HPP

#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cmath>
#include <sstream>
#include <fstream>

using namespace std;


struct rec_time {
    long int date;
    long int seconds;
    double milli_seconds;
};


//convert string to rec_time
rec_time StringToTime(string str) {
    rec_time result;
    string tmp;
    string::iterator it;
    
    //get data from string
    for (it = str.begin(); (*it) != ':'; ++it) {
        tmp += *it;
    }
    long int date = stoi(tmp.c_str());
    
    //get seconds from string
    tmp.clear();
    for (; (*it) != '.'; ++it) {
        if( (*it)!=':' ) tmp +=*it;
    }
    long int seconds = stoi(tmp.c_str());
    
    //get milli seconds from string
    tmp.clear();
    for (; it != str.end(); ++it) {
        tmp += *it;
    }
    double m_seconds = stof(tmp.c_str());
    
    result.date = date;
    result.seconds = seconds;
    result.milli_seconds = m_seconds;
    
    return result;
}


struct record {
    rec_time time;
    float price;
    int volume;
};


//convert record string to record vector
vector<record> StringToRecord (string source){
    vector<record> result_vec;
    string buffer;
    string::iterator it = source.begin();
    
    //skip the first line for the possible broken data during slicing the big data file
    while (it != source.end() && (*it) != '\n') {
        ++it;
    }
    ++it;
    
    for (; it != source.end(); ++it) {
        
        try {
            //get the first line
            while ((*it) != '\n') {
                buffer += *it;
                ++it;
            }
            
            string::iterator it1 = buffer.begin();
            string tmp_str;
            record tmp_record;
            
            //get record time
            while ((*it1) != ',') {
                tmp_str += *it1;
                ++it1;
            }
            tmp_record.time = StringToTime(tmp_str);
            tmp_str.clear();
            
            //get price
            ++it1;
            while ((*it1) != ',') {
                tmp_str += *it1;
                ++it1;
            }
            tmp_record.price =stof(tmp_str);
            tmp_str.clear();
            
            //get volume
            ++it1;
            while (it1 != buffer.end()) {
                tmp_str += *it1;
                ++it1;
            }
            tmp_record.volume = stoi(tmp_str);
            
            result_vec.push_back(tmp_record);
            buffer.clear();
            
            if(it == source.end()) return result_vec;
            
        } catch (...) {
            cout << "This line is wrong format! Ignor it!" << endl;
        }
        
    }
    
    return result_vec;
}

string RecordToString(vector<record> & record_source) {
    stringstream ss;
    long int tmp1, tmp2, tmp3, tmp4, tmp;
    
    for (auto it = record_source.begin(); it != record_source.end(); ++it) {
        tmp = (*it).time.seconds;
        tmp1 = tmp / 10000;
        tmp %= 10000;
        tmp2 = tmp / 100;
        tmp3 = tmp % 100;
        tmp4 = (long int)((*it).time.milli_seconds * 1000000);
        
        ss << (*it).time.date << ':' << tmp1 << ':' << tmp2 << ':' << tmp3 << '.' << tmp4 << ',' << (*it).price << ',' << (*it).volume << '\n';
    }
    
    return ss.str();
}


//Test if record 2 is later than record 1
bool Later (record & r1, record & r2){
    if(r1.time.date != r2.time.date) {
        return (r1.time.date < r2.time.date) ? true : false;
    }
    else if (r2.time.seconds != r2.time.seconds) {
        return (r1.time.seconds < r2.time.seconds) ? true : false;
    }
    else {
        return (r1.time.milli_seconds < r2.time.milli_seconds) ? true : false;
    }
}


//sort the slide window
void SortWindow(vector<record> & rec_vector, int start, int end) {
    int window_size = end - start + 1;
    for (int i = 0; i < window_size; ++i)
    {
        for (int j = i+1; j< window_size; ++j)
            if (Later(rec_vector[j], rec_vector[i]))
            {
                swap(rec_vector[j], rec_vector[i]);
            }
    }
}


//scrub the date to two parts
pair<vector<record>, vector<record> > ScrubRecord(vector<record> & source, int window_size, double tolerance = 0.30) {
    vector<record> signal;
    vector<record> noise;
    vector<record> a = source;
    int n = a.size();

    if (window_size > n) {
        window_size = n - 4;
    }
    
    //initialize slide window
    SortWindow(source, 0, window_size - 1);
    
    double sum = 0.0;
    double squared_sum = 0.0;
    double mean = 0.0;
    double std_err = 0.0;
    
    //sort sliding window
    for (int i = window_size; i < n; ++i) {
        
        if (Later(a[i-1], source[i])) {
        }
        
        else {
            SortWindow(a, i - window_size, i);
        }
    }
    
    for (int i=1; i<n-1; ++i) {
        if ( (a[i].price < min(a[i-1].price, a[i+1].price) / 10.0) || (a[i].price > max(a[i-1].price, a[i+1].price) * 10.0)) {
            a[i].price = (a[i-1].price + a[i+1].price ) / 2.0;
        }
        
        sum += a[i].price;
        squared_sum +=(a[i].price - mean) * (a[i].price - mean);
//        cout << "i: " << i << ", price:" << source[i].price << ", squared_sum:" << squared_sum << endl;
        
    }
    
    mean = sum / (n-2.0);
//    cout << "mean: " << mean << endl;
    
    std_err = sqrt(squared_sum / (n-2.0)) ;
//    cout << "std_err: " << std_err << endl;
    
    //seperate signal and noise
    for (int i =0; i < n; ++i) {
        if (abs(source[i].price - mean) > tolerance * std_err) {
            noise.push_back(source[i]);
        }
        else signal.push_back(source[i]);
    }
    
    return make_pair(signal, noise);
}




bool NormalTest(vector<record> & source) {
    long N = source.size();
    vector<double> return_vec(N-1, 0);
    
    //initialize return vector
    for (int i=0; i<N; ++i) {
        return_vec[i] = (source[i+1].price - source[i].price) / source[i].price;
    }
    
    vector<long double> moment(4,0);    //first to forth moments
    vector<long double> sum(4,0);       //sum, sum_square, sum_cubic, sum_quad
    
    for (int i=0; i < N-1; ++i) {
        sum[0] += return_vec[i];
        sum[1] += pow(return_vec[i], 2.0);
        sum[2] += pow(return_vec[i], 3.0);
        sum[3] += pow(return_vec[i], 4.0);
    }
    
    moment[0] = sum[0] / (double)(N-1.0);
    moment[1] = sum[1] - 2.0 * moment[0] * sum[0] - (N-1.0) * pow(moment[0], 2.0);
    moment[2] = sum[2] - 3.0 * sum[1] * moment[0] + 3.0 * pow(moment[0], 2.0) * sum[0] - pow(moment[0], 3.0);
    moment[3] = sum[3] - 4.0 * moment[1] * sum[2] + 6.0 * pow(moment[0], 2.0) * sum[1] -4.0 * pow(moment[0], 3.0) * sum[0] + (N-1.0) * pow(moment[0], 4.0);
    
    //JB test, return ture if pass normality test， confident level 99%
//    cout << (N-1.0)*pow(moment[2], 2.0)/6.0 + (N-1.0)*pow(moment[3]-3.0, 2.0)/24.0 <<endl;
    return ( (N-1.0)*pow(moment[2], 2.0)/6.0 + (N-1.0)*pow(moment[3]-3.0, 2.0)/24.0 ) < 9.21;
}



#endif /* scrub_HPP */
