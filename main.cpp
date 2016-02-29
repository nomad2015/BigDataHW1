//
//  main.cpp
//  9898BDiFHW1
//
//  Created by Xiaobo on 2/20/16.
//  Copyright © 2016 Xiaobo. All rights reserved.
//

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
#include "mpi.h"

#include "scrub.hpp"

using namespace std;

int main(int argc, const char * argv[]) {
    
//    ifstream in_file;
//    in_file.open("/Users/xiaobohe/Desktop/c++/9898_BDiF_HW/HW1/9898BDiFHW1/9898BDiFHW1/data10k.txt");
//    
//    string str, tmp;
//    int counter = 0;
//    
//    while (!in_file.eof()) {
//        getline(in_file, tmp);
//        tmp += '\n';
//        str += tmp;
//        ++counter;
//    }
//
//    cout << counter << endl;
//   
//    vector<record> record_test = StringToRecord(str);
//
//  
//    ofstream signal_file, noise_file;
//    signal_file.open("/Users/xiaobohe/Desktop/c++/9898_BDiF_HW/HW1/9898BDiFHW1/9898BDiFHW1/signal.txt");
//    noise_file.open("/Users/xiaobohe/Desktop/c++/9898_BDiF_HW/HW1/9898BDiFHW1/9898BDiFHW1/noise.txt");
//    
//    vector<record> signal_test, noise_test;
//    
//    pair<vector<record>, vector<record>> result_test;
//    result_test = ScrubRecord(record_test, 50);
//    
//    signal_test = result_test.first;
//    noise_test = result_test.second;
//
//    string signal_str = RecordToString(signal_test);
//    string noise_str = RecordToString(noise_test);
//    
//    signal_file << signal_str;
//    noise_file << noise_str;
//    
//    signal_file.close();
//    noise_file.close();
//
//    
//    cout << "size signal: " << signal_test.size() << ", size noise: " << noise_test.size() <<endl;
//    
//    
//    //test normality
//    bool normality = NormalTest(signal_test);
//    cout << normality << endl;
    
    
    
    
    
    

    
    int rank, nodes;
    
    //start MPI
    MPI_File fh;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nodes);
    MPI_File_open(MPI_COMM_WORLD,argv[1],MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
    
    MPI_Offset FILESIZE;//=stoi(argv[1]);
    MPI_File_get_size(fh, &FILESIZE);
    
    int size_buf = FILESIZE/nodes;
    char buff[size_buf/sizeof(char)];
    
    //read into buff
    clock_t t_scrub_start = clock();    //time to start reading
    MPI_File_read_at(fh, rank*size_buf, buff, size_buf/sizeof(char), MPI_BYTE, &status);
    
    //convert string to record vector
    vector<record> record_vec = StringToRecord(buff);
    MPI_File_close(&fh);
    
    //scrub data into signal and noise
    vector<record> signal_vec, noise_vec;
    pair<vector<record>, vector<record>> result_vec;
    result_vec = ScrubRecord(record_test, 50);
    signal_vec = result_vec.first;
    noise_vec = result_vec.second;
    
    //time to finish scrubing
    clock_t t_scrub_end = clock();
    
    //output signal into signal.txt file
    string signal_str = RecordToString(signal_vec);
    MPI_Offset output_offset = signal_str.size();
    MPI_Offset * send_offset = new long long;   //
    *send_offset=output_offset;
    long long * rbuf = (long long *)malloc(size*sizeof(long long));
    MPI_Allgather( send_offset, 1, MPI_LONG, rbuf, 1, MPI_LONG, MPI_COMM_WORLD);
    
    MPI_Offset cumulative_offset = 0;
    for (int i=0; i < rank; ++i) {
        cumulative_offset += rbuf[i];
    }
    MPI_File fh_out;
    MPI_File_open(MPI_COMM_WORLD, "signal.txt", MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_out);
    char * p_signal = new char[signal_str.size()];
    strcpy(p_signal,signal_str.c_str());
    
    MPI_File_write_at(fh_out, cumulative_offset, p_signal, output_offset, MPI_BYTE, &status);
    MPI_File_close(&fh_out);
    
    //output noise into noise.txt file
    string noise_str = RecordToString(noise_vec);
    output_offset=noise_str.size();
    delete send_offset;
    send_offset = new long long;
    *send_offset=output_offset;
    free(rbuf);
    rbuf = (long long *)malloc(size*sizeof(long long));
    MPI_Allgather( send_offset, 1, MPI_LONG, rbuf, 1, MPI_LONG, MPI_COMM_WORLD);
    cumulative_offset = 0;
    for (int i = 0; i < rank; ++i){
        cumulative_offset += rbuf[i];
    }
    MPI_File_open(MPI_COMM_WORLD, "noise.txt", MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_out);
    char * p_noise=new char[noise_str.size()];
    strcpy(p_noise,noise_str.c_str());
    MPI_File_write_at(fh_out, cumulative_offset, p_noise, output_offset, MPI_BYTE, &status);
    MPI_File_close(&fh_out);
    
    
    
    //Test normality with one node
    //time to start normality test
    clock_t t_nor_start = clock();
    if (rank == 0) {
        
        bool normality_result = NormalTest(signal_vec);
        
        stringstream ss;
        if (normality_result) {
            ss << "Normality test passed with JB method!" <<endl;
        }
        else ss << "Normality test did not pass with JB method!" << endl;
        
        ofstream out_file;
        out_file.open("NormalityTest.txt");
        out_file << ss.str();
        out_file.close();
    }
    //time to finish normality test
    clock_t t_nor_end = clock();
    
    
    //write log file
    if (rank == 0) {
        ofstream log_file;
        log_file.open("log.txt");
        stringstream log_ss;
        log_ss << "At node 0, time to read data to scrub: " << t_scrub_start << endl << "Time to finish scrubing: " << t_scrub_end << endl << "It took " << t_scrub_end - t_scrub_start << " to scrub the data." << endl << endl << "Time to start normality test: " << t_nor_start <<endl << "Time to finish nornality test: " << t_nor_end << endl << "It took " << t_nor_end - t_nor_start << " to do normality test."  << endl;
        log_file << log_ss.str();
        log_file.close();
    }
    
    MPI_Finalize();
    return 0;

}
