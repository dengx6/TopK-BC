 /*
 * topkbc.h
 *
 *  Created on: Apr 3, 2023
 *      Author: dengxin
 */

#ifndef TOPKBC_H_
#define TOPKBC_H_

#include "../edge/edge.h"
#include "../gstream/gstream.h"
#include "../utils/utils.h"

//#define DEBUG_TRACK
//#define TOPK_RESULT_OUTPUT_FILE
//#define INITIAL_TIME_OUTPUT_FILE
#define GLOBAL_PRUNE
#define LOCAL_PRUNE
#define GLOBAL_COMMENT 1

//#define SORT_BY_DEGREE
//#define SORT_BY_CORE
#define SORT_BY_QWEIGHTS
#define SORT_BY_QWEIGHTS_DY

#define SAMPLE_RATE 0.01
#define SZ(x) ((int)x.size())
#define max(x,y) ((x)>(y)?(x):(y))
#define min(x,y) ((x)<(y)?(x):(y))

typedef long long lint;

#ifdef SORT_BY_DEGREE
class VertexDegree{
public:
    unsigned vertex;
    unsigned degree;
    inline bool operator < (const VertexDegree &other) const {
        return degree > other.degree || (degree == other.degree && vertex < other.vertex);
    }

    VertexDegree();
    VertexDegree(unsigned v, unsigned d);
    ~VertexDegree();
};
#endif


class VertexQWeight{
public:
	unsigned vertex;
	Wtype qweight;
	inline bool operator < (const VertexQWeight &other) const{
		return qweight > other.qweight || (qweight == other.qweight && vertex < other.vertex);
	}
	
	VertexQWeight();
	VertexQWeight(unsigned vertex, Wtype qweight);
	~VertexQWeight();
	
};

class Neighbor{
public:
	unsigned op_vertex;
	Edge* related_edge;

	inline bool operator < (const Neighbor &other) const {
		return op_vertex < other.op_vertex;
	}

	inline bool operator == (const Neighbor &other) const {
		return op_vertex == other.op_vertex;
	}

	Neighbor();
	Neighbor(unsigned _op_vertex, Edge* _related_edge);
	~Neighbor();
};

class SmallBclique{
public:
	unsigned R_vertex;
	Wtype Lweight;
	lint min_ts;

	inline bool operator < (const SmallBclique &other) const {
		return R_vertex < other.R_vertex;
	}

	SmallBclique();
	SmallBclique(unsigned _R_vertex, Wtype _Lweight, lint min_ts);
	~SmallBclique();

};

class pqbclique{
public:
	Edge* related_edge;
	lint min_timestamp;
	vector<Edge*> edges_of_bclique;

	pqbclique();
	pqbclique(Edge* related_edge, lint min_timestamp, vector<Edge*> &_edges_of_bclique);
	~pqbclique();
};

class topkbc{
public:
	topkbc(unsigned _p, unsigned _q, unsigned _winsz, unsigned _stride, string _data_file_path, unsigned long long _update_times, 
			unsigned _k, unsigned _sort_strategy, unsigned _average_time_span, unsigned long long _initial_time,
			unsigned long long _final_time);
	~topkbc();
	void run();
	void insert(Edge* _newedge);
	void remove();
	void construct_global_index(Edge* _newedge, bool insert_or_delete);
	bool construct_local_graph(unsigned l, vector<SmallBclique> &_R, vector<unsigned> &_CL);
	bool construct_local_graph(Edge* _newedge, bool insert_or_delete);
	bool trim_graph_by_core(bool insert_or_delete);
	double estimate_cost(unsigned side);
	void sort_vertices(unsigned strategy);
	void collect_two_hop_adj();
	void serach_Topk_bcliques();
	void PQbclique(unsigned l);
	void find_combinations(vector< vector<Neighbor> > &temp_R, vector<SmallBclique> &temp_aux_R, vector< vector<Neighbor> > &temp_R_comb,
			vector<SmallBclique> &temp_aux_R_comb, SmallBclique &op_vert_sbclique,vector<Neighbor> &op_vert_Neighbors,
			unsigned beg_offset, unsigned curr_depth, unsigned total_depth);
	Wtype replace_Topk_and_update_densityUbMap();
	void remove_expired_result();
	void search_alternate_biclique();

public:
	lint tmp_ts;
	lint num_edges;
	lint updated_edges_insert;
	lint updated_edges_delete;
	unsigned long long update_count;
	bool is_same_time;
	bool insert_or_delete; 
	
	unsigned long long win_start_time;
	unsigned long long win_end_time;

	unsigned sliding_counter;
	Wtype cur_maximum_density_ub;
	Wtype next_maximum_density_ub;
	unsigned n_vertices_local[2];
	unsigned num_vertices;
	unsigned t;
	unsigned op_t;
	unsigned vert;
	unsigned op_vert;
	Wtype Minimum_density_needed;

	// time
	double time_output, time_initial, time_monitor, time_insert, time_delete;
	double tbegin_monitor, tend_monitor;
	double t_begin,t_end;
	double tbegin_initial, tend_initial;
	double tbegin_insert, tend_insert;
	double tbegin_delete, tend_delete;
	lint average_space;
	lint total_space;
	
	Edge* cur_edge;

	vector<unsigned> tmp_vec;

	map<unsigned, vector<Neighbor>> global_one_hop_neighbor[2];
	
	map<unsigned, vector<Neighbor>> local_one_hop_neighbor[2];//sorted by Vertex lexicographical orde
	
	map<unsigned, Edge*> cur_edge_one_hop_neighbor[2];

	VertexQWeight cur_edge_qweight[2];

	map<unsigned, vector<unsigned>> two_hop_adj_map;
	
	vector<unsigned> sorted_vertices;
	
	vector<VertexQWeight> verQWeights;

	unordered_map<unsigned, unsigned> vertices_rank;
	
#ifndef SORT_BY_QWEIGHTS
	unordered_map<unsigned, Wtype> vertices_qweights;
#endif
	
	vector< vector<unsigned> > L;

	vector< vector< vector<Neighbor> > > R;

	vector< vector<SmallBclique>> aux_vec_R;

	vector< vector<unsigned> > Cand_of_L;

	vector< map<unsigned, vector<Neighbor>>> candidate_one_hop_neighbor;

	vector< vector<VertexQWeight>> candidate_verQWeights;

	vector< map<unsigned, vector<unsigned>>> candidate_two_hop_adj;

	vector< map<unsigned, unsigned>> candidate_vertices_rank;

	map<Edge*, Wtype> map1_of_densityUb;

	multimap<Wtype, Edge*, greater<Wtype>> map2_of_densityUb;

	multimap<Wtype, pqbclique, greater<Wtype>> map_of_Topk_res;

	vector<Edge*> snapshot;

public:
	//util
	void print_local_one_hop_neighbor();
	void print_sorted_vertices();
	void print_local_two_hop_neighbor();
	void print_cur_edge_one_hop_neighbor();
	void print_Topk_weight_bicliques();
	void print_map_of_edgeRank();
	void print_edges();
	long getMemoryUse();

	string make_outfile_name(string dataname, unsigned pv, unsigned qv, lint update_num, unsigned clique_num, unsigned edge_num);
	static vector<unsigned> intersection(vector<unsigned> &vec1, vector<unsigned> &vec2, int offset1, int offset2);
    vector<Neighbor> intersection(bool insert_or_delete, vector<Neighbor> &vec1, vector<Neighbor> &vec2, int offset1, int offset2);
    vector<vector<Neighbor>> intersection(unsigned l, vector<vector<Neighbor>> &vec1, vector<Neighbor> &vec2, vector<SmallBclique> &vec3, int offset1, int offset2);
	static vector<Neighbor> intersection(vector<SmallBclique> &vec1, vector<Neighbor> &vec2, map<unsigned, vector<unsigned>> &R_one_hop_neig, unsigned vertex, int offset1, int offset2);
	static vector<Neighbor> merge(vector<Neighbor> &vec1, vector<Neighbor> &vec2);
	static unsigned intersection_count(vector<unsigned> &vec1, vector<unsigned> &vec2);

private:
	unsigned p;
	unsigned q;
	unsigned k;
	unsigned winsz;
	unsigned stride;
	string data_file_path;
	unsigned long long update_times;
	bool anchor_left;
	unsigned sort_strategy;
	unsigned average_time_span;
	unsigned long long initial_time;
	unsigned long long final_time;
};

#endif /* TOPKBC_H_ */