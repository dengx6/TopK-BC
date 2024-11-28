#ifndef _GSTREAM_H_
#define _GSTREAM_H_

#include "../edge/edge.h"
class gstream{
public:
	gstream(string _data_path);
	~gstream();
	bool hasnext();
	Edge * next();
	int size();

private:
	ifstream fin;
	long long edge_number;
};


#endif
