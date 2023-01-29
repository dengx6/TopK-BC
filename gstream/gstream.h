#ifndef _GSTREAM_H_
#define _GSTREAM_H_

#include "edge.h"
class gstream{
public:
	gstream(string _data_path);
	~gstream();
	bool hasnext();
	Edge * next();


private:
	ifstream fin;

};


#endif
