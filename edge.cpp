#include "edge.h"

Edge::Edge(int _u, int _v, int _weight, long long _timestamp):u(_v),v(_v),weight(_weight),timestamp(_timestamp)
{
		
}

Edge::Edge(string _str_e){
	stringstream _ss(_str_e);
	_ss >> this->u >> this->v;
	_ss >> this->weight >> this->timestamp;
}
