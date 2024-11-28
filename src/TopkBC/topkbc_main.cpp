//============================================================================
// Name        : topkbc-main.cpp
// Author      : dengxin
// Version     : 2024.03.20
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "topkbc.h"
#define LOCAL_PRUNE

int main(int argc, char* argv[]) {
	int data_set = 0;
	string data_home = "xxx";
	string str_dat[1] = {
		"xxx.txt" // dataset name
	};
	string _dataset = data_home + str_dat[data_set];
	
	double begin_clock, end_clock, elapsed_time;
	unsigned _p = 4;
	unsigned _q = 4;
	unsigned _window = 10000;//time window
	unsigned _update_times = 10000;
	unsigned _average_time_span = 1;
	unsigned _stride = 1;
	unsigned _k = 32;
	unsigned sort_strategy = 2;
	unsigned long long _initial_time = 1;
	unsigned long long _final_time = 1000000;
	
	cout << "dataset = " << _dataset << endl;
	cout << "p = " << _p << endl;
	cout << "q = " << _q << endl;
	cout << "window = " << _window << endl;
	cout << "update_times= " << _update_times  <<endl;
	cout << "stride = " << _stride << endl;
	cout << "k = " << _k << endl;

	topkbc tbc(_p, _q, _window, _stride, _dataset, _update_times, _k, sort_strategy, _average_time_span, _initial_time, _final_time);
	begin_clock = clock();
	tbc.run();
	end_clock = clock();
	elapsed_time = (end_clock - begin_clock) / CLOCKS_PER_SEC;
	
	double time_update = tbc.time_insert + tbc.time_delete;
	long long total_edge_update = tbc.updated_edges_insert + tbc.updated_edges_delete;
	
	cout << "" << endl;
	cout << "dataset = " << _dataset << endl;
	cout << "p = " << _p << "q = " << _q << "k = " << _k << endl;
	cout << "window = " << _window << endl;
	cout << "Our final total runing time: using " << elapsed_time << " seconds." << endl;
	cout << "Throughput = " << (double)(total_edge_update/time_update) << " (edges/sec)." << endl;
	return 0;
}
