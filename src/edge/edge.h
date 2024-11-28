#ifndef _EDGE_H_
#define _EDGE_H_

#include <iostream>
#include <set>
#include <vector>
#include <fstream>
#include <sstream>
#include <queue>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <algorithm>
#include <iterator>
#include <string.h>
//#include <pthread.h>
#include <unistd.h>
#include <sys/time.h>
#include <cstdio>
#include <string>
#include <cstring>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <time.h>
#include <cmath>
#include <random>
#include <sys/resource.h>
using namespace std;

typedef double Wtype;

class Edge{
public:
	Edge(unsigned _u, unsigned _v, Wtype _weight, long long _timestamp, long long number);
	Edge(string _str_e, long long number);
	unsigned u;
	unsigned v;
	Wtype weight;
	long long timestamp;
	long long number;

	inline bool operator < (const Edge &other) const{
		return (u < other.u) || (u == other.u && v < other.v);
	}

	inline bool operator == (const Edge &other) const{
		return u == other.u && v == other.v;
	}
};
#endif
