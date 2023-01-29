//============================================================================
// Name        : main.cpp
// Author      : dengxin
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================
#include "edge.h"
#include "basesolution.h"
/* 
 * argv0 : exe
 * argv1 : dataset
 * argv2 : p
 * argv3 : q
 * argv4 : winsz
 * 
 * */
int main(int argc, char* argv[])
{
/*
#ifdef DEBUG_TRACK
	
#endif

#ifdef THREAD_LOG
#endif
	if(argc < 5){
		cout << "err argc" << endl;
		exit(0);
	}
	cout << "num of argc: " << argc << endl;

	long long _update_times = 5000;
	int k = 10;
	string _dataset;
	_dataset = string(argv[1]);
	
	int _p, _q, _window;
	{
		stringstream _ss;
		for(int i = 2; i < argc; i ++) 
			_ss << argv[i] << " ";
		_ss >> _p >> _q>> _window;
	}
	cout << "dataset=" << _dataset << endl;
	cout << "p=" << _p << endl;
	cout << "q=" << _q << endl;
	cout << "window=" << _window << endl;
*/
	double begin_clock, end_clock, elapsed_time;
	int _p = 6;
	int _q = 4;
	int _window = 100;
	string _dataset = "D:\\Programma project\\eclipseWorkspace\\top-k_test\\datasets\\test\\rel_sort.rating";
	int _update_times = 1000;
	int k = 2;
#ifdef DEBUG
	cout << "dataset = " << _dataset << endl;
	cout << "p = " << _p << endl;
	cout << "q = " << _q << endl;
	cout << "window = " << _window << endl;
	cout << "_update_times= " << _update_times  <<endl;
	cout << "k = " << k << endl;
	getchar();
#endif
	basesolution bsolut(_p, _q, _window, _dataset, _update_times, k);
	begin_clock = clock();
	bsolut.run();
	end_clock = clock();
	elapsed_time = (end_clock - begin_clock) / CLOCKS_PER_SEC;
	cout << "total runing time: using " << elapsed_time << " seconds." << endl;
	cout << "The top-k result weight as follow :" << endl;
	for(auto itr : bsolut.Resultmap){
		cout << itr.first << endl;
	}
	return 0;
}
