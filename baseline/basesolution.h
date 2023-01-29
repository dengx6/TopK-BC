#ifndef _BASESOLUTION_H_
#define _BASESOLUTION_H_

#define DEBUG 1
#define GLOBAL_COMMENT 1

#define MAX_P 20

#include "edge.h"
#include "../gstream/gstream.h"
//顶点邻居
struct neighbor{
	int x;
	int weight;
	long long timestamp;
	
	neighbor(int _x, int _weight, long long _timestamp):x(_x), weight(_weight), timestamp(_timestamp){};
	~neighbor(){};
};

class basesolution{
public:
	basesolution(int _p, int _q, int _winsz , string _data_path, long long _update_times, int k);
	~basesolution();
	void run();
	void update(Edge* _newedge);
	void insert(Edge* _newedge);
	void remove();
	void collect_2hop_neighbors(Edge* _newedge);
	void init_partial_biclique(Edge* _newedge);
	void neighbor_intersection(int l, int &_c, vector<int>* _R);//update R
	void update_Candidate(int l, int i, vector<int>* _R, vector<int>* _C);//update C
	void recursionTree_construction(int l, vector<int> &_L);
	void find_combinations(vector< vector<int> > &combs, int cur_picked, int min_subscript, vector<int> &_R, vector<int> &comb);
	//update resultmap and edgemap;
	void rm_old_result(long long _timestamp);//remove expired result from resultmap and edgemap
	void add_new_result(vector<int> &_L, vector<int>* _R); //add new result to resultmap and edgemap

	int snapshot_size();
	
	vector<Edge*>::iterator tmp_itr;//迭代器
	long long tmp_ts;
	bool ts_is_same;
	vector<Edge*> snapshot;

	set<int> neighbor_of_u;
	set<int> neighbor_of_v;
	multimap<int,set<Edge*>, greater<int>> Resultmap;//top-k结果集，来源于辅助结构，key是权重之和W(B),降序排序
	multimap<int,set<Edge*>, greater<int>> Edgemap;//辅助结构，用来存储当前所有的(p,q)-biclique，删除了top-k (p,q)-biclique
	
	map<int, int> commom_neighbor_cnt;//记录插入边的顶点u,v与二跳邻居之间共同邻居数量,目前只记录顶点u的二跳邻居
	map<int,set<Edge*>> one_hop_tb;//u或v所有二跳邻居的邻居表,目前只记录顶点u的二跳邻居
	//map<int,set<int>> two_hop_tb;//顶点u或v的二跳邻居表

	vector<int> L;//L of partial biclique
	vector<int>* r[MAX_P];
	vector<int>* c[MAX_P];
private:
	int p;
	int q;
	int k;
	int winsz;
	string data_path;
	long long update_times;
	long long update_cnt;
};



#endif
