#include "edge.h"

Edge::Edge(unsigned _u,  unsigned _v, Wtype _weight, long long _timestamp, long long _number)
{
		this->u = _u;
		this->v = _v;
		this->weight = _weight;
		this->timestamp = _timestamp;
		this->number = _number;
}

Edge::Edge(string _str_e, long long _number){
	stringstream _ss(_str_e);
	_ss >> this->u >> this->v;
	_ss >> this->weight >> this->timestamp;
	this->number = _number;
}
