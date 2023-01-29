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
#include <unistd.h>
#include <sys/time.h>

using namespace std;
class Edge{
public:
	Edge(int _u, int _v, int _weight, long long _timestamp);
	Edge(string _str_e);
	int u;
	int v;
	int weight;
	long long timestamp;

	bool operator < (const Edge &other) const{
		return (u < other.u) || (u == other.u && v < other.v);
	}
};
#endif
