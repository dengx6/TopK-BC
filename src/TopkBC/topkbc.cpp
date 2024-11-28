#include "topkbc.h"

#ifdef SORT_BY_DEGREE
VertexDegree::VertexDegree(): vertex(-1), degree(-1){}

VertexDegree::VertexDegree(unsigned v, unsigned d): vertex(v), degree(d){}

VertexDegree::~VertexDegree(){}
#endif

bool cmpByWeightDec(Neighbor &v1, Neighbor &v2){
    return v1.related_edge->weight > v2.related_edge->weight ||
		   (v1.related_edge->weight == v2.related_edge->weight && v1.op_vertex < v2.op_vertex);
}

bool cmpByLWeightsDec(SmallBclique &v1, SmallBclique &v2){

    return v1.Lweight > v2.Lweight ||
		   (v1.Lweight == v2.Lweight && v1.R_vertex < v2.R_vertex);
}

bool cmpByPairSecond(pair<unsigned, Edge*> &v1, pair<unsigned, Edge*> &v2){
	return v1.second->weight > v2.second->weight || (v1.second->weight == v2.second->weight && v1.first < v2.first);
}

VertexQWeight::VertexQWeight(): vertex(0), qweight(0){}

VertexQWeight::VertexQWeight(unsigned _vertex, Wtype _qweight): vertex(_vertex), qweight(_qweight){}

VertexQWeight::~VertexQWeight(){}

Neighbor::Neighbor(): op_vertex(0), related_edge(nullptr){}

Neighbor::Neighbor(unsigned _op_vertex, Edge* _related_edge){
	op_vertex = _op_vertex;
	related_edge = _related_edge;
}

Neighbor::~Neighbor(){}

SmallBclique::SmallBclique(): R_vertex(0), Lweight(0), min_ts(0){};

SmallBclique::SmallBclique(unsigned _R_vertex, Wtype _Lweight, lint _min_ts){
	R_vertex = _R_vertex;
	Lweight = _Lweight;
	min_ts = _min_ts;
}
SmallBclique::~SmallBclique(){};

pqbclique::pqbclique(): related_edge(nullptr),min_timestamp(0){}

pqbclique::pqbclique(Edge* _related_edge, lint _min_timestamp, vector<Edge*> &_edges_of_bclique){
	related_edge = _related_edge;
	min_timestamp = _min_timestamp;
	edges_of_bclique = _edges_of_bclique;
}

pqbclique::~pqbclique(){}

topkbc::topkbc(unsigned _p, unsigned _q, unsigned _winsz, unsigned _stride, string _data_file_path, unsigned long long _update_times, unsigned _k, 
		unsigned _sort_strategy, unsigned _average_time_span, unsigned long long _initial_time, unsigned long long _final_time)
{
		this->p = _p;
		this->q = _q;
		this->winsz = _winsz;
		this->data_file_path = _data_file_path;
		this->update_times =_update_times;
		this->stride = _stride;
		this->update_count = 0;
		this->k = _k;
		this->sort_strategy = _sort_strategy;
		
		this->average_time_span = _average_time_span;
		this->initial_time = _initial_time;
		this->final_time = _final_time;
		this->average_space = 0;
		this->total_space = 0;
		
		this->anchor_left = 0;
		this->tmp_ts = 0;
		this->is_same_time = 0;
		this->num_edges = 0;
		this->sliding_counter = 0;
		this->time_output = 0;
		this->time_initial = 0;
		this->time_monitor = 0;
		this->time_insert = 0;
		this->time_delete = 0;
		
		this->snapshot.clear();
		this->cur_edge = nullptr;
}

topkbc::~topkbc(){
	for(int i = 0; i < (int)this->snapshot.size();i ++){
		delete this->snapshot[i];
		this->snapshot[i] = NULL;//先释放指针，后置NULL
	}
	this->snapshot.clear();
}


void topkbc::run(){	
	this->tbegin_monitor = clock();
	this->tbegin_initial = clock();
	
	this->win_start_time = this->initial_time;
	this->win_end_time = this->initial_time + this->winsz * this->average_time_span;
	this->updated_edges_insert = 0;
	this->updated_edges_delete = 0;
    double time_update = 0.0;
    long long update_edges_total = 0;
    
	gstream gs(this->data_file_path);
    Edge* newedge = nullptr;
    while(gs.hasnext()){
    	newedge = gs.next();
    	if(newedge->timestamp < this->win_end_time){
    		if(this->update_count >= 1){
    			this->tbegin_insert = clock();
    			this->updated_edges_insert ++;
    		}
    		
    		this->insert(newedge);
    		if(this->update_count >= 1){
    			this->tend_insert = clock();
    			this->time_insert += (this->tend_insert - this->tbegin_insert) / CLOCKS_PER_SEC;
    		}
    	
    	}
    	else{
    		if(this->update_count >= 1){
    			this->tbegin_delete = clock();
    		}
    		this->remove();
    		if(this->update_count >= 1){
    			this->tend_delete = clock();
    			this->time_delete += (this->tend_delete - this->tbegin_delete) / CLOCKS_PER_SEC;
    		}
    		
    		this->update_count ++;
    		if(this->update_count == 1){
    			// compute initial time
    			this->tend_initial = clock();
    			this->time_initial = (this->tend_initial - this->tbegin_initial) / CLOCKS_PER_SEC;
    			cout << "|------------Initialization completed---------------|" << endl;
    		}
    		this->tend_monitor = clock();
    		this->time_monitor = (this->tend_monitor - this->tbegin_monitor) / CLOCKS_PER_SEC;//compute current running time
    		time_update = this->time_insert + this->time_delete;
    		update_edges_total = this->updated_edges_insert + this->updated_edges_delete;
    		
    		if(this->update_count <= this->update_times){
    			this->total_space += mem::getValue();
    			this->average_space = this->total_space / this->update_count;
    		}

    		this->t_begin = clock();
    		this->print_Topk_weight_bicliques();
    					
			this->t_end = clock();
			this->time_output += (this->t_end - this->t_begin) / CLOCKS_PER_SEC;
    					
			// update time interval, the window will slide
    		this->win_end_time = this->win_end_time + this->average_time_span;// interval close
    		this->win_start_time = this->win_start_time + this->average_time_span; // interval open
    		if(this->win_end_time >= this->final_time){
    			return;
    		}
    		
			while(newedge->timestamp >= this->win_end_time){
        		if(this->update_count >= 1){
        			this->tbegin_delete = clock();
        		}
    			this->remove();
        		if(this->update_count >= 1){
        			this->tend_delete = clock();
        			this->time_delete += (this->tend_delete - this->tbegin_delete) / CLOCKS_PER_SEC;
        		}
        		this->update_count ++;
        		
        		this->tend_monitor = clock();
        		this->time_monitor = (this->tend_monitor - this->tbegin_monitor) / CLOCKS_PER_SEC;//compute current running time
        		time_update = this->time_insert + this->time_delete;
        		update_edges_total = this->updated_edges_insert + this->updated_edges_delete;
        		
        		if(this->update_count <= this->update_times){
        			this->total_space += mem::getValue();
        			this->average_space = this->total_space / this->update_count;
        		}
        		
        		this->t_begin = clock();
    			this->print_Topk_weight_bicliques();
    			    			
    			this->t_end = clock();
    			this->time_output += (this->t_end - this->t_begin) / CLOCKS_PER_SEC;
    			
    			// update time interval, the window will slide
        		this->win_end_time = this->win_end_time + this->average_time_span;// interval close
        		this->win_start_time = this->win_start_time + this->average_time_span; // interval open
        		if(this->win_end_time >= this->final_time){
        			return;
        		}
			}
			
    		if(this->update_count >= 1){
    			this->tbegin_insert = clock();
    		}
			this->insert(newedge);
			this->updated_edges_insert ++;
    		if(this->update_count >= 1){
    			this->tend_insert = clock();
    			this->time_insert += (this->tend_insert - this->tbegin_insert) / CLOCKS_PER_SEC;
    		}
    	}
    }
}

void topkbc::insert(Edge* _newedge){
	this->snapshot.push_back(_newedge);
	this->num_edges ++;
	this->construct_global_index(_newedge, 1);
	this->cur_edge = _newedge; 
	this->insert_or_delete = true;
	
	if(this->construct_local_graph(_newedge, 1)){
		this->serach_Topk_bcliques();
	}
	//cout << "insert finished" <<endl;
}

void topkbc::remove(){
	unsigned old_topk_size = this->map_of_Topk_res.size();
	this->remove_expired_result();
	if(old_topk_size >= k && this->map_of_Topk_res.size() < this->k){
		this->search_alternate_biclique();
	}
	//cout << "remove finished" <<endl;
}

void topkbc::construct_global_index(Edge* _newedge, bool insert_or_delete){
	unsigned u = _newedge->u;
	unsigned v = _newedge->v;

	if(insert_or_delete){
		Neighbor u_neighbor(v, _newedge);
		Neighbor v_neighbor(u, _newedge);
		this->global_one_hop_neighbor[0][u].push_back(u_neighbor);
		this->global_one_hop_neighbor[1][v].push_back(v_neighbor);
	}
	else{
		auto u_neighbor_itr = this->global_one_hop_neighbor[0][u].begin();
		auto v_neighbor_itr = this->global_one_hop_neighbor[1][v].begin();
		for(unsigned i = 0; i < this->global_one_hop_neighbor[0][u].size(); i++){
			if(this->global_one_hop_neighbor[0][u][i].op_vertex == v){
				this->global_one_hop_neighbor[0][u].erase(u_neighbor_itr + i);
				break;
			}
		}
		for(unsigned i = 0; i < this->global_one_hop_neighbor[1][v].size(); i++){
			if(this->global_one_hop_neighbor[1][v][i].op_vertex == u){
				this->global_one_hop_neighbor[1][v].erase(v_neighbor_itr + i);
				break;
			}
		}
	}

}

bool topkbc::construct_local_graph(unsigned l, vector<SmallBclique> &_R, vector<unsigned> &_CL){
	//cout << "IN construct_local_graph_of_serach_Topk " << endl;
	//step 1: compute candidate one hopneighbor and use (q, p)-core to graph reduction

	this->candidate_vertices_rank[l-1].clear();
	this->candidate_one_hop_neighbor[l-1].clear();
	map<unsigned, vector<unsigned>> R_one_hop_neighbor;
	//cout << "_CL size = " << _CL.size() << endl;
	for(unsigned i = 0; i < _CL.size(); i++){
		unsigned vertex = _CL[i];
		vector<Neighbor> neig_vec = this->local_one_hop_neighbor[1 - this->anchor_left][vertex];
		vector<Neighbor> com_neig_vec = this->intersection(_R, neig_vec, R_one_hop_neighbor, vertex, 0, 0);
		if(com_neig_vec.size() >= this->op_t){
			//cout << "vertex = " << vertex <<endl;
			this->candidate_one_hop_neighbor[l-1].insert({vertex, com_neig_vec});
#ifndef SORT_BY_QWEIGHTS
		this->candidate_vertices_rank[l-1].insert({vertex, i});
#endif
		}
	}

	//step 2: vertices sort by Qweight again
	unsigned num_cand = this->candidate_one_hop_neighbor[l-1].size();
	if(num_cand != _CL.size()){
		cerr << "err in construct local graph" << endl;
		getchar();
	}
	if(l + num_cand < this->t - 1){
		return false;
	}
	//this->candidate_vertices_rank[l-1].clear();
	this->candidate_verQWeights[l-1].resize(num_cand);
	Wtype *QWeights = new Wtype[num_cand]();
	unsigned j = 0;

	for(auto itr : this->candidate_one_hop_neighbor[l-1]){
		sort(itr.second.begin(), itr.second.end(), cmpByWeightDec);
		for(unsigned i = 0, m = 0; i < (this->op_t - 1) && (m < itr.second.size()); m++){
			unsigned neighbor = itr.second[m].op_vertex;
			if(neighbor != this->op_vert){
				QWeights[j] += itr.second[m].related_edge->weight;
				i++;
			}
		}
		//cout << "this->cur_edge_one_hop_neighbor[this->anchor_left][itr.first]->weight = " << this->cur_edge_one_hop_neighbor[this->anchor_left][itr.first]->weight <<endl;
		QWeights[j] += this->cur_edge_one_hop_neighbor[this->anchor_left][itr.first]->weight;
		VertexQWeight vq(itr.first, QWeights[j]);
		this->candidate_verQWeights[l-1][j] = vq;
		j++;
	}
	sort(this->candidate_verQWeights[l-1].begin(), this->candidate_verQWeights[l-1].end());
#ifdef SORT_BY_QWEIGHTS
	this->Cand_of_L[l-1].resize(num_cand);
	for(unsigned i = 0; i < num_cand; i++){
		//cout << "this->candidate_verQWeights[l-1][i].vertex = " << this->candidate_verQWeights[l-1][i].vertex <<endl;
		this->Cand_of_L[l-1][i] = this->candidate_verQWeights[l-1][i].vertex;
		this->candidate_vertices_rank[l-1].insert({this->candidate_verQWeights[l-1][i].vertex, i});
	}
#endif
	delete[] QWeights;

	//step 3: dynantic compute two hop neighbor
#ifndef SORT_BY_QWEIGHTS
	this->candidate_two_hop_adj[l-1].clear();
	map<unsigned, vector<unsigned> > aux_array_two_neig;
	for(unsigned i = 0; i < num_cand; i++){
		unsigned vertex = _CL[i];
		for(unsigned j = 0; j < this->candidate_one_hop_neighbor[l-1][vertex].size(); j++){
			unsigned neig = this->candidate_one_hop_neighbor[l-1][vertex][j].op_vertex;
			if(neig == this->op_vert){continue;}
			for(unsigned k = 0; k < R_one_hop_neighbor[neig].size(); k++){
				unsigned two_hop_neighbor = R_one_hop_neighbor[neig][k];
				auto rank_itr = this->candidate_vertices_rank[l-1].find(two_hop_neighbor);
				if(rank_itr != this->candidate_vertices_rank[l-1].end()){
					unsigned two_hop_neig_rank = rank_itr->second;
					if(two_hop_neig_rank > i){
						aux_array_two_neig[two_hop_neighbor].push_back(neig);
					}
				}
			}
		}
		for(auto itr : aux_array_two_neig){
			unsigned two_hop_neighbor = itr.first;
			unsigned num_common_neig = itr.second.size();
			if(num_common_neig >= this->op_t - 1){
				this->candidate_two_hop_adj[l-1][vertex].push_back(two_hop_neighbor);
			}
		}
		aux_array_two_neig.clear();
	}
#endif

#ifdef SORT_BY_QWEIGHTS
	this->candidate_two_hop_adj[l-1].clear();
	map<unsigned, vector<unsigned> > aux_array_two_neig;
	for(unsigned i = 0; i < num_cand; i++){
		unsigned vertex = this->candidate_verQWeights[l-1][i].vertex;
		for(unsigned j = 0; j < this->candidate_one_hop_neighbor[l-1][vertex].size(); j++){
			unsigned neig = this->candidate_one_hop_neighbor[l-1][vertex][j].op_vertex;
			if(neig == this->op_vert){continue;}
			for(unsigned k = 0; k < R_one_hop_neighbor[neig].size(); k++){
				unsigned two_hop_neighbor = R_one_hop_neighbor[neig][k];
				auto rank_itr = this->candidate_vertices_rank[l-1].find(two_hop_neighbor);
				if(rank_itr != this->candidate_vertices_rank[l-1].end()){
					unsigned two_hop_neig_rank = rank_itr->second;
					if(two_hop_neig_rank > i){
						aux_array_two_neig[two_hop_neighbor].push_back(neig);
					}
				}
			}
		}
		//compute num of commom neighbor
		for(auto itr : aux_array_two_neig){
			unsigned two_hop_neighbor = itr.first;
			unsigned num_common_neig = itr.second.size();
			if(num_common_neig >= this->op_t - 1){
				//unsigned two_hop_neig_rank = this->candidate_vertices_rank[l-1][two_hop_neighbor];
				this->candidate_two_hop_adj[l-1][vertex].push_back(two_hop_neighbor);
			}
		}
//		sort(this->candidate_two_hop_adj[l-1][vertex].begin(), this->candidate_two_hop_adj[l-1][vertex].end());//ascend sort
//		for(unsigned i = 0; i < this->candidate_two_hop_adj[l-1][vertex].size(); i++){
//			unsigned rank = this->candidate_two_hop_adj[l-1][vertex][i];
//			this->candidate_two_hop_adj[l-1][vertex][i] = this->candidate_verQWeights[l-1][rank].vertex;
//		}
		aux_array_two_neig.clear();
	}
#endif
	//cout << "OUT construct_local_graph_of_serach_Topk" <<endl;
	return true;
}

bool topkbc::construct_local_graph(Edge* _newedge, bool insert_or_delete){//可增加判断，局部子图是否有效
	if(!this->trim_graph_by_core(insert_or_delete)){
		//cout << "No results because the graph is pruned by core" << endl;
		return false;
	}

	unsigned u = _newedge->u;
	unsigned v = _newedge->v;
	this->n_vertices_local[0] = this->local_one_hop_neighbor[0].size();
	this->n_vertices_local[1] = this->local_one_hop_neighbor[1].size();
	this->num_vertices = this->n_vertices_local[0] + this->n_vertices_local[1];

#ifdef DEBUG_TRACK
	//print_local_one_hop_neighbor();
#endif

	if(n_vertices_local[0] < this->p - 1 || n_vertices_local[1] < this->q - 1){
		//cout << "No results because the graph is pruned by core" << endl;
		return false;
	}


	//step 2:compute qweights of currently updating edge
	vector<pair<unsigned, Edge*>> vec1(this->cur_edge_one_hop_neighbor[0].begin(), this->cur_edge_one_hop_neighbor[0].end());
	vector<pair<unsigned, Edge*>> vec2(this->cur_edge_one_hop_neighbor[1].begin(), this->cur_edge_one_hop_neighbor[1].end());

	sort(vec1.begin(), vec1.end(), cmpByPairSecond);
	sort(vec2.begin(), vec2.end(), cmpByPairSecond);
	Wtype qweight[2];
	memset(qweight, 0, sizeof(qweight));
	for(unsigned i = 0, j = 0; i < this->n_vertices_local[1] && j < this->q-1; i++){
		unsigned neighbor = vec1[i].first;
		if(neighbor != v){
			qweight[0] += vec1[i].second->weight;
			j++;
		}
	}
	qweight[0] += this->cur_edge->weight;
	for(unsigned i = 0, j = 0; i < this->n_vertices_local[0] && j < this->p-1; i++){
		unsigned neighbor = vec2[i].first;
		if(neighbor != u){
			qweight[1] += vec2[i].second->weight;
			j++;
		}
	}
	qweight[1] += this->cur_edge->weight;
	this->cur_edge_qweight[0].vertex = u;
	this->cur_edge_qweight[0].qweight = qweight[0];
	this->cur_edge_qweight[1].vertex = v;
	this->cur_edge_qweight[1].qweight = qweight[1];

	//compute estimate_cost of local graph
	if(estimate_cost(0) < estimate_cost(1)){
		this->anchor_left = true;
	}
	else{
		this->anchor_left = false;
	}
	//cout << "anchor_left = " << this->anchor_left << endl;

	this->t = this->anchor_left == 1 ? this->p : this->q;
	this->op_t = this->anchor_left == 1 ? this->q : this->p;
	this->vert = this->anchor_left == 1 ? this->cur_edge->u : this->cur_edge->v;
	this->op_vert = this->anchor_left == 1 ? this->cur_edge->v : this->cur_edge->u;

	this->sort_vertices(this->sort_strategy);
#ifdef DEBUG_TRACK
	//this->print_sorted_vertices();
#endif
	//cout << "finish sorting vertices" << endl;

	this->collect_two_hop_adj();

#ifdef DEBUG_TRACK
	//this->print_local_two_hop_neighbor();
#endif
	return true;
}

bool topkbc::trim_graph_by_core(bool insert_or_delete){
	/*
	 * 注意插入与删除时构建局部子图的方式不同
	 * 删除时,诱发的局部子图每条边的edge number均小于搜索边的edge number
	 */
	unsigned u = this->cur_edge->u;
	unsigned v = this->cur_edge->v;
	vector<Neighbor> intersect_ans;
	unsigned num_u_neig = this->global_one_hop_neighbor[0][u].size();
	unsigned num_v_neig = this->global_one_hop_neighbor[1][v].size();
	if(num_u_neig  < this->q || num_v_neig < this->p){
		//cout << "currently updated edge degree less than q or p, do not need to use (q,p)-core" <<endl;
		return false;
	}

	this->cur_edge_one_hop_neighbor[0].clear();
	this->cur_edge_one_hop_neighbor[1].clear();
	this->local_one_hop_neighbor[0].clear();
	this->local_one_hop_neighbor[1].clear();
	sort(this->global_one_hop_neighbor[0][u].begin(), this->global_one_hop_neighbor[0][u].end());
	sort(this->global_one_hop_neighbor[1][v].begin(), this->global_one_hop_neighbor[1][v].end());

	unsigned start_idx[2], end_idx[2];
	unsigned **to_remove_vertices = new unsigned*[2];
	memset(start_idx, 0 , sizeof(start_idx));
	memset(end_idx, 0 , sizeof(end_idx));
	to_remove_vertices[0] = new unsigned[num_v_neig]();
	to_remove_vertices[1] = new unsigned[num_u_neig]();
	set<unsigned> removed[2];

	//step 1: collect the initial vertices with degree less than q in left and p in right
	if(insert_or_delete){
		//cout << "num_u_neig = " << num_u_neig <<endl;
		//cout << "num_v_neig = " << num_v_neig <<endl;
		for(unsigned i = 0, j = 0; i < num_u_neig || j < num_v_neig; i++, j++){
			//cout << "i = " << i << " j = " << j <<endl;
			//cout << "this->global_one_hop_neighbor[0][u][1].op_vertex = " << this->global_one_hop_neighbor[0][u][1].op_vertex <<endl;
			if(i < num_u_neig){
				unsigned u_neighbor = this->global_one_hop_neighbor[0][u][i].op_vertex;
				//cout << "u_neighbor = " << u_neighbor << endl;
				if(u_neighbor != v){
					sort(this->global_one_hop_neighbor[1][u_neighbor].begin(), this->global_one_hop_neighbor[1][u_neighbor].end());
					intersect_ans = this->intersection(1, this->global_one_hop_neighbor[1][u_neighbor], this->global_one_hop_neighbor[1][v], 0 , 0);
					this->local_one_hop_neighbor[1].insert({u_neighbor, intersect_ans});
					if(intersect_ans.size() < this->p){
						//cout << "to remove u_neighbor = " << u_neighbor <<endl;
						to_remove_vertices[1][end_idx[1]++] = u_neighbor;
						removed[1].insert(u_neighbor);
					}
					else{
						this->cur_edge_one_hop_neighbor[0].insert({u_neighbor, this->global_one_hop_neighbor[0][u][i].related_edge});
					}
				}
			}

			if(j < num_v_neig){
				unsigned v_neighbor = this->global_one_hop_neighbor[1][v][j].op_vertex;
				//cout << "v_neighbor = " << v_neighbor << endl;
				if(v_neighbor != u){
					sort(this->global_one_hop_neighbor[0][v_neighbor].begin(), this->global_one_hop_neighbor[0][v_neighbor].end());
					intersect_ans = this->intersection(1, this->global_one_hop_neighbor[0][v_neighbor], this->global_one_hop_neighbor[0][u], 0 , 0);
					this->local_one_hop_neighbor[0].insert({v_neighbor, intersect_ans});
					if(intersect_ans.size() < this->q){
						//cout << "to remove v_neighbor = " << v_neighbor <<endl;
						to_remove_vertices[0][end_idx[0]++] = v_neighbor;
						removed[0].insert(v_neighbor);
					}
					else{
						this->cur_edge_one_hop_neighbor[1].insert({v_neighbor, this->global_one_hop_neighbor[1][v][j].related_edge});
					}
				}
			}
		}
	}
	else{
		//cout << "num_u_neig = " << num_u_neig <<endl;
		//cout << "num_v_neig = " << num_v_neig <<endl;
		for(unsigned i = 0, j = 0; i < num_u_neig || j < num_v_neig; i++, j++){
			//cout << "i = " << i << " j = " << j <<endl;
			if(i < num_u_neig){
				unsigned u_neighbor = this->global_one_hop_neighbor[0][u][i].op_vertex;
				lint edge_number = this->global_one_hop_neighbor[0][u][i].related_edge->number;
				if(u_neighbor != v && edge_number < this->cur_edge->number){
					sort(this->global_one_hop_neighbor[1][u_neighbor].begin(), this->global_one_hop_neighbor[1][u_neighbor].end());
					intersect_ans = this->intersection(0, this->global_one_hop_neighbor[1][u_neighbor], this->global_one_hop_neighbor[1][v], 0 , 0);
					this->local_one_hop_neighbor[1].insert({u_neighbor, intersect_ans});
					if(intersect_ans.size() < this->p){
						//cout << "to remove u_neighbor = " << u_neighbor <<endl;
						to_remove_vertices[1][end_idx[1]++] = u_neighbor;
						removed[1].insert(u_neighbor);
					}
					else{
						this->cur_edge_one_hop_neighbor[0].insert({u_neighbor, this->global_one_hop_neighbor[0][u][i].related_edge});
					}

				}
			}
			//cout << "this->global_one_hop_neighbor[1][v][0].op_vertex = "  << this->global_one_hop_neighbor[1][v][0].op_vertex;
			if(j < num_v_neig){
				unsigned v_neighbor = this->global_one_hop_neighbor[1][v][j].op_vertex;
				lint edge_number = this->global_one_hop_neighbor[1][v][j].related_edge->number;
				//cout << "v_neighbor = " << v_neighbor <<endl;
				//cout << "edge_number = " << edge_number <<endl;
				//cout << "this->cur_edge->number = " << this->cur_edge->number <<endl;
				if(v_neighbor != u && edge_number < this->cur_edge->number){
					sort(this->global_one_hop_neighbor[0][v_neighbor].begin(), this->global_one_hop_neighbor[0][v_neighbor].end());
					intersect_ans = this->intersection(0, this->global_one_hop_neighbor[0][v_neighbor], this->global_one_hop_neighbor[0][u], 0 , 0);
					this->local_one_hop_neighbor[0].insert({v_neighbor, intersect_ans});
					//cout << "2-this->local_one_hop_neighbor[0] size = " << this->local_one_hop_neighbor[0].size() << endl;
					if(intersect_ans.size() < this->q){
						to_remove_vertices[0][end_idx[0]++] = v_neighbor;
						removed[0].insert(v_neighbor);
					}
					else{
						this->cur_edge_one_hop_neighbor[1].insert({v_neighbor, this->global_one_hop_neighbor[1][v][j].related_edge});
					}
				}
			}
		}
	}
	if(removed[0].size() == 0 && removed[1].size() == 0){
		for(unsigned i = 0; i < 2; i++){
			delete[] to_remove_vertices[i];
		}
		delete[] to_remove_vertices;
		return true;
	}

	//cout << "0-this->local_one_hop_neighbor[1] size = " << this->local_one_hop_neighbor[1].size() << endl;
	unsigned initial_removed_u = end_idx[0];
	unsigned initial_removed_v = end_idx[1];

	//step 2: recursively remove all vertices with degree less than q in left and p in right, i.e., (q,p)-core
	while(start_idx[0] != end_idx[0] || start_idx[1] != end_idx[1]){
		//cout << "end_idx[0] = " << end_idx[0] << endl;
		if(start_idx[0] != end_idx[0]){
			unsigned u_vertex = to_remove_vertices[0][start_idx[0]++];
			//cout << "remove u_vertex = " << u_vertex << endl;
			for(unsigned i = 0; i < this->local_one_hop_neighbor[0][u_vertex].size(); i++){
				unsigned other = this->local_one_hop_neighbor[0][u_vertex][i].op_vertex;
				if(other == v){
					if(start_idx[0] <= initial_removed_u){
						continue;
					}
					if(this->cur_edge_one_hop_neighbor[1].find(u_vertex) != this->cur_edge_one_hop_neighbor[1].end()){
						this->cur_edge_one_hop_neighbor[1].erase(u_vertex);
					}
					if(this->cur_edge_one_hop_neighbor[1].size() < this->p - 1){
						for(unsigned i = 0; i < 2; i++){
							delete[] to_remove_vertices[i];
						}
						delete[] to_remove_vertices;
						return false;
					}
				}
				else{
					if(removed[1].find(other) == removed[1].end()){
						Neighbor neig(u_vertex, nullptr);
						auto del_vertex = lower_bound(this->local_one_hop_neighbor[1][other].begin(), this->local_one_hop_neighbor[1][other].end(), neig);
						if(del_vertex != this->local_one_hop_neighbor[1][other].end() && del_vertex->op_vertex == u_vertex){
							this->local_one_hop_neighbor[1][other].erase(del_vertex);
						}
						if(this->local_one_hop_neighbor[1][other].size() < this->p){
							to_remove_vertices[1][end_idx[1]++] = other;
							removed[1].insert(other);
						}
					}
				}
			}
		}
		//cout << "1-this->local_one_hop_neighbor[1] size = " << this->local_one_hop_neighbor[1].size() << endl;
		if(start_idx[1] != end_idx[1]){
			//cout << "end_idx[1] = " << end_idx[1] << endl;
			unsigned v_vertex = to_remove_vertices[1][start_idx[1]++];
			//cout << "remove v_vertex = " << v_vertex << endl;
			for(unsigned i = 0; i < this->local_one_hop_neighbor[1][v_vertex].size(); i++){
				unsigned other = this->local_one_hop_neighbor[1][v_vertex][i].op_vertex;
				//cout << "other = " << other << endl;
				if(other == u){
					if(start_idx[1] <= initial_removed_v){
						continue;
					}
					if(this->cur_edge_one_hop_neighbor[0].find(v_vertex) != this->cur_edge_one_hop_neighbor[0].end()){
						this->cur_edge_one_hop_neighbor[0].erase(v_vertex);
					}
					else{
						cout << "cur_edge_one_hop_neighbor[0]" << endl;
						getchar();
					}
					if(this->cur_edge_one_hop_neighbor[0].size() < this->q - 1){
						for(unsigned i = 0; i < 2; i++){
							delete[] to_remove_vertices[i];
						}
						delete[] to_remove_vertices;
						return false;
					}
				}
				else{
					if(removed[0].find(other) == removed[0].end()){
						Neighbor neig(v_vertex, nullptr);
						auto del_vertex = lower_bound(this->local_one_hop_neighbor[0][other].begin(), this->local_one_hop_neighbor[0][other].end(), neig);
						if(del_vertex != this->local_one_hop_neighbor[0][other].end() && del_vertex->op_vertex == v_vertex){
							this->local_one_hop_neighbor[0][other].erase(del_vertex);
						}
						else{
							cout << "error!!! not find del_vertex" << endl;
							getchar();
						}

						if(this->local_one_hop_neighbor[0][other].size() < this->q){
							to_remove_vertices[0][end_idx[0]++] = other;
							removed[0].insert(other);
						}
					}
				}
			}
		}
	}

	//step 3: update local one hop neighbors(此处可优化)
	auto itr1 = this->local_one_hop_neighbor[0].begin();
	auto itr2 = this->local_one_hop_neighbor[1].begin();
	//cout << "2-this->local_one_hop_neighbor[1] size = " << this->local_one_hop_neighbor[1].size() << endl;
	for(; itr1 != this->local_one_hop_neighbor[0].end() || itr2 != this->local_one_hop_neighbor[1].end();){
		if(itr1 != this->local_one_hop_neighbor[0].end()){
			if(removed[0].find(itr1->first) != removed[0].end()){
				this->local_one_hop_neighbor[0].erase(itr1++);
			}
			else{
				itr1++;
			}
		}

		if(itr2 != this->local_one_hop_neighbor[1].end()){
			if(removed[1].find(itr2->first) != removed[1].end()){
				this->local_one_hop_neighbor[1].erase(itr2++);
			}
			else{
				itr2++;
			}
		}

	}
	//cout << "3-this->local_one_hop_neighbor[1] size = " << this->local_one_hop_neighbor[1].size() << endl;
	for(unsigned i = 0; i < 2; i++){
		delete[] to_remove_vertices[i];
	}
	delete[] to_remove_vertices;

	//cout << "after trimming by core" << endl;
	return true;
}

double topkbc::estimate_cost(unsigned side){
	vector<pair<unsigned, vector<Neighbor>>> vec(this->local_one_hop_neighbor[side].begin(), this->local_one_hop_neighbor[side].end());
	unsigned vertex = side == 0 ? this->cur_edge->u : this->cur_edge->v;
	unsigned op_vertex = side == 0 ? this->cur_edge->v : this->cur_edge->u;
	map<unsigned, unsigned> aux_array_two_neig;
	srand (time(NULL));

	unsigned num_rounds = ceil(num_vertices * SAMPLE_RATE);

	lint total_two_hop_deg = 0;
	lint max_two_hop_deg = 0;
	
	//unsigned *common_neig_map = new unsigned[this->num_vertices]();
	//unsigned *aux_array_two_neig = new unsigned[this->num_vertices]();
	//unsigned offset = side == 0 ? 0 : n_vertices[0];
	unsigned common_neig_threshold = side == 0 ? this->q - 1 : this->p - 1;

	for(unsigned r = 0; r < num_rounds; r++){
		lint estimated_two_hop_deg = 0;
		//unsigned i = rand() % n_vertices[side] + offset;
		unsigned i = rand() % n_vertices_local[side];
		unsigned i_vertex = vec[i].first;

		for(unsigned j = 0; j < vec[i].second.size(); j++){
			unsigned neighbor = vec[i].second[j].op_vertex;
			if(neighbor == op_vertex){continue;}
			for(unsigned k = 0; k < this->local_one_hop_neighbor[1-side][neighbor].size(); k++){
				unsigned two_hop_neighbor = this->local_one_hop_neighbor[1-side][neighbor][k].op_vertex;
				if(two_hop_neighbor == vertex){continue;}
				if(two_hop_neighbor < i_vertex){
					if(aux_array_two_neig.find(two_hop_neighbor) != aux_array_two_neig.end()){
						aux_array_two_neig[two_hop_neighbor]++;
					}
					else{
						aux_array_two_neig[two_hop_neighbor] = 1;
					}
				}
				else{
					break;
				}
			}
		}

		for(auto itr : aux_array_two_neig){
			if(itr.second >= common_neig_threshold){
				estimated_two_hop_deg += 1;
			}
		}
		aux_array_two_neig.clear();

		max_two_hop_deg = max_two_hop_deg < estimated_two_hop_deg ? estimated_two_hop_deg : max_two_hop_deg;
		total_two_hop_deg = total_two_hop_deg + estimated_two_hop_deg * this->n_vertices_local[side];
	}

	total_two_hop_deg = total_two_hop_deg / num_rounds;
	//cout << "estimated_total_two_hop_deg[" << side << "]=" << total_two_hop_deg << ", estimated_avg_two_hop_deg[" << side << "]=" << total_two_hop_deg / this->n_vertices_local[side] << ", max_two_hop_deg[" << side << "]=" << max_two_hop_deg << endl;

	lint avg_two_hop_deg = max(2, total_two_hop_deg / this->n_vertices_local[side]);  // we let avg_two_hop_deg be at least 2 to deal with corner case
	lint pq_value = side == 0? this->p - 1: this->q - 1;
	double totalCost = total_two_hop_deg * pow(avg_two_hop_deg, pq_value-2);
	//cout << "totalCost = " << totalCost << endl;

	return totalCost;
}

void topkbc::sort_vertices(unsigned strategy){
	 switch (strategy) {
		case 0: {//degree vertex order
			//cout << "degree vertex order is applied..." << endl;
#ifdef SORT_BY_DEGREE
			this->sorted_vertices.clear();
			this->vertices_rank.clear();
			this->vertices_qweights.clear();
			this->verQWeights.clear();
			this->verQWeights.resize(this->n_vertices_local[1- this->anchor_left]);

			vector<VertexDegree> verdegs(this->n_vertices_local[1- this->anchor_left]);
			Wtype *QWeights = new Wtype[this->n_vertices_local[1- this->anchor_left]]();
			unsigned j = 0;
			//compute the degree and qweights of vertices
			for(auto itr : this->local_one_hop_neighbor[1 - this->anchor_left]){
				sort(itr.second.begin(), itr.second.end(), cmpByWeightDec);
				for(unsigned i = 0, m = 0; i < (this->op_t - 1) && (m < itr.second.size()); m++){
					unsigned neighbor = itr.second[m].op_vertex;
					if(neighbor != this->op_vert){
						QWeights[j] += itr.second[m].related_edge->weight;
						i++;
					}
				}
				QWeights[j] += this->cur_edge_one_hop_neighbor[this->anchor_left][itr.first]->weight;
				this->vertices_qweights.insert({itr.first, QWeights[j]});

				VertexDegree vd(itr.first, itr.second.size());
				verdegs[j] = vd;
				j++;
			}
			sort(verdegs.begin(), verdegs.end());
			this->sorted_vertices.resize(this->n_vertices_local[1- this->anchor_left]);
			for(unsigned i = 0; i < this->n_vertices_local[1- this->anchor_left]; i++){
				VertexQWeight vq(verdegs[i].vertex, this->vertices_qweights[verdegs[i].vertex]);
				this->sorted_vertices[i] = verdegs[i].vertex;
				this->verQWeights[i] = vq;
				this->vertices_rank[verdegs[i].vertex] = i;

			}
			delete[] QWeights;
#endif
			break;
		}
		case 1: {// core vertex order
			//cout << "core vertex order is applied..." << endl;
#ifdef SORT_BY_CORE
			this->sorted_vertices.clear();
			this->vertices_rank.clear();
			this->vertices_qweights.clear();
			this->verQWeights.clear();
			this->verQWeights.resize(this->n_vertices_local[1- this->anchor_left]);
			this->sorted_vertices.resize(this->n_vertices_local[1- this->anchor_left]);

			map<unsigned, unsigned> vertices_deg[2];
			Wtype *QWeights = new Wtype[this->n_vertices_local[1- this->anchor_left]]();
			unsigned j = 0;

			//step 1: collect degree and qweights of vertices
			for(auto itr : this->local_one_hop_neighbor[this->anchor_left]){
				unsigned vertex = itr.first;
				unsigned degree = itr.second.size();
				vertices_deg[1].insert({vertex, degree});
			}
			for(auto itr : this->local_one_hop_neighbor[1 - this->anchor_left]){
				unsigned vertex = itr.first;
				unsigned degree = itr.second.size();
				vertices_deg[0].insert({vertex, degree});
				sort(itr.second.begin(), itr.second.end(), cmpByWeightDec);
				for(unsigned i = 0, m = 0; i < (this->op_t - 1) && (m < degree); m++){
					unsigned neighbor = itr.second[m].op_vertex;
					if(neighbor != this->op_vert){
						QWeights[j] += itr.second[m].related_edge->weight;
						i++;
					}
				}
				QWeights[j] += this->cur_edge_one_hop_neighbor[this->anchor_left][itr.first]->weight;
				this->vertices_qweights.insert({itr.first, QWeights[j]});
				j++;
			}

			//step 2: compute core vertex order
			unsigned start_idx[2], end_idx[2];
			unsigned **to_remove_vertices = new unsigned*[2];
			memset(start_idx, 0 , sizeof(start_idx));
			memset(end_idx, 0 , sizeof(end_idx));
			to_remove_vertices[0] = new unsigned[this->n_vertices_local[1- this->anchor_left]]();
			to_remove_vertices[1] = new unsigned[this->n_vertices_local[this->anchor_left]]();
			set<unsigned> removed[2];

			unsigned idx = this->n_vertices_local[1- this->anchor_left];
			unsigned p_value = this->t;
			unsigned q_value = this->op_t;
			while(idx > 0){
				for(unsigned i = 0; i < end_idx[0]; i++){
					to_remove_vertices[0][i] = 0;
				}
				for(unsigned j = 0; j < end_idx[1]; j++){
					to_remove_vertices[1][j] = 0;
				}
				start_idx[0] = end_idx[0] = 0;
				start_idx[1] = end_idx[1] = 0;

				for(auto itr : vertices_deg[0]){
					unsigned vertex = itr.first;
					unsigned degree = itr.second;
					auto set_itr = removed[0].find(vertex);
					if(degree < q_value && set_itr == removed[0].end()){
						to_remove_vertices[0][end_idx[0]++] = vertex;
						removed[0].insert(vertex);
						this->sorted_vertices[--idx] = vertex;
						this->vertices_rank.insert({vertex, idx});
						VertexQWeight vq(vertex, this->vertices_qweights[vertex]);
						this->verQWeights[idx] = vq;
					}
				}
				for(auto itr : vertices_deg[1]){
					unsigned vertex = itr.first;
					unsigned degree = itr.second;
					auto set_itr = removed[1].find(vertex);
					if(degree < p_value && set_itr == removed[1].end()){
						to_remove_vertices[1][end_idx[1]++] = vertex;
						removed[1].insert(vertex);
					}
				}
				 while(start_idx[0] != end_idx[0] || start_idx[1] != end_idx[1]){
					 if(start_idx[0] != end_idx[0]){
						 unsigned vertex = to_remove_vertices[0][start_idx[0]++];

						 for(int i = 0; i < vertices_deg[0][vertex]; i++){
							 unsigned other = this->local_one_hop_neighbor[1-this->anchor_left][vertex][i].op_vertex;
							 if(other == this->op_vert){continue;}
							 for(int j = 0; j < vertices_deg[1][other]; j++){
								 unsigned other_neig = this->local_one_hop_neighbor[this->anchor_left][other][j].op_vertex;
								 if(other_neig == vertex){
									 vertices_deg[1][other]--;
									 break;
								 }
							 }

							 auto set_itr = removed[1].find(other);
							 if(vertices_deg[1][other] < p_value && set_itr == removed[1].end()){
								 to_remove_vertices[1][end_idx[1]++] = other;
								 removed[1].insert(other);
							 }
						 }
						 vertices_deg[0][vertex] = 0;
						 continue;
					 }

					 if(start_idx[1] != end_idx[1]){
						 unsigned vertex = to_remove_vertices[1][start_idx[1]++];

						 for(int i = 0; i < vertices_deg[1][vertex]; i++){
							 unsigned other = this->local_one_hop_neighbor[this->anchor_left][vertex][i].op_vertex;
							 if(other == this->vert){continue;}
							 for(int j = 0; j < vertices_deg[0][other]; j++){
								 unsigned other_neig = this->local_one_hop_neighbor[1-this->anchor_left][other][j].op_vertex;
								 if(other_neig == vertex){
									 vertices_deg[0][other]--;
									 break;
								 }
							 }
							 auto set_itr = removed[0].find(other);
							 if(vertices_deg[0][other] < q_value && set_itr == removed[0].end()){
								 to_remove_vertices[0][end_idx[0]++] = other;
								 removed[0].insert(other);
								 this->sorted_vertices[--idx] = other;
								 this->vertices_rank.insert({other, idx});
								 VertexQWeight vq(other, this->vertices_qweights[other]);
								 this->verQWeights[idx] = vq;
							 }
						 }
						 vertices_deg[1][vertex] = 0;
					 }
				 }
				 p_value++;
				 q_value++;
			}

			for(unsigned i = 0; i < 2; i++){
				delete[] to_remove_vertices[i];
			}
			delete[] to_remove_vertices;
			delete[] QWeights;
#endif
			break;
		}
		case 2: { // qweights vertex order
			//cout << "qweights vertex order is applied..." << endl;
#ifdef SORT_BY_QWEIGHTS
			this->sorted_vertices.clear();
			this->vertices_rank.clear();
			this->verQWeights.clear();
			this->verQWeights.resize(this->n_vertices_local[1- this->anchor_left]);

			Wtype *QWeights = new Wtype[this->n_vertices_local[1- this->anchor_left]]();//debug
			unsigned j = 0;

			for(auto itr : this->local_one_hop_neighbor[1- this->anchor_left]){
				sort(itr.second.begin(), itr.second.end(), cmpByWeightDec);
				for(unsigned i = 0, m = 0; i < (this->op_t - 1) && (m < itr.second.size()); m++){
					unsigned neighbor = itr.second[m].op_vertex;
					if(neighbor != this->op_vert){
						QWeights[j] += itr.second[m].related_edge->weight;
						i++;
					}
				}
				QWeights[j] += this->cur_edge_one_hop_neighbor[this->anchor_left][itr.first]->weight;
				VertexQWeight vq(itr.first, QWeights[j]);
				this->verQWeights[j] = vq;
				j++;
			}
			sort(this->verQWeights.begin(), this->verQWeights.end());
			for(unsigned i = 0; i < verQWeights.size(); i++){
				this->sorted_vertices.push_back(verQWeights[i].vertex);
				this->vertices_rank[verQWeights[i].vertex] = i;
			}

			delete[] QWeights;
#endif
			break;
		}
		default:{
			cout << "Please select the following vertex ordering: 0-random, 1-degree, 2-core" << endl;
			return;
		}
	}
}

void topkbc::collect_two_hop_adj(){
	this->two_hop_adj_map.clear();
	map<unsigned, vector<unsigned> > aux_array_two_neig;

	for(unsigned i = 0; i < this->sorted_vertices.size(); i++){
		unsigned vertex = this->sorted_vertices[i];
		for(unsigned j = 0; j < this->local_one_hop_neighbor[1 - this->anchor_left][vertex].size(); j++){
			unsigned neighbor = this->local_one_hop_neighbor[1 - this->anchor_left][vertex][j].op_vertex;
			if(neighbor == this->op_vert){continue;}
			for(unsigned k = 0; k < this->local_one_hop_neighbor[this->anchor_left][neighbor].size(); k++){
				unsigned two_hop_neighbor = this->local_one_hop_neighbor[this->anchor_left][neighbor][k].op_vertex;
				if(two_hop_neighbor == this->vert){continue;}
				if(this->vertices_rank[two_hop_neighbor] > i){//Rank at the following vertex
					aux_array_two_neig[two_hop_neighbor].push_back(neighbor);
				}
			}
		}
		for(auto itr : aux_array_two_neig){
			unsigned two_hop_neighbor = itr.first;
			unsigned num_common_neig = itr.second.size();
			if(num_common_neig >= this->op_t - 1){
				this->two_hop_adj_map[vertex].push_back(two_hop_neighbor);
			}
		}
		aux_array_two_neig.clear();
	}
}

void topkbc::serach_Topk_bcliques(){
	this->cur_maximum_density_ub = 0;
	if(this->insert_or_delete){
		if(this->map_of_Topk_res.size() < this->k){
			this->Minimum_density_needed = 0;
			//cout << "insert: initial Minimum density needed = " << this->Minimum_density_needed << endl;
		}
		else{
			this->Minimum_density_needed = this->map_of_Topk_res.rbegin()->first;
			//cout << "insert: initial Minimum density needed = " << this->Minimum_density_needed << endl;
		}
	}
	else{
		if(this->map_of_Topk_res.size() < this->k){
			this->Minimum_density_needed = this->next_maximum_density_ub;
		}
		else{
			this->Minimum_density_needed = this->map_of_Topk_res.rbegin()->first;
		}
		//cout << "delete: initial Minimum density needed = " << this->Minimum_density_needed << endl;
	}

	//start to serach
	for(unsigned i = 0; i < this->n_vertices_local[1 - this->anchor_left]; i++){
		/**
		 * global pruning
		 *
		 */
		// cout << "start to global pruning" <<endl;
		if(i == this->n_vertices_local[1 - this->anchor_left] - this->t + 2){
			break;
		}
		Wtype global_density_ub = 0;
		unsigned m = 0;

#ifndef SORT_BY_QWEIGHTS
		vector<VertexQWeight> ver_qweight;
		ver_qweight.assign(this->verQWeights.begin() + i + 1, this->verQWeights.end());
		
		sort(ver_qweight.begin(), ver_qweight.end());
		for(unsigned j = 0; m < this->t - 2 && j < ver_qweight.size(); j++,m++){
			global_density_ub += ver_qweight[m].qweight;
		}
		if(m == this->t - 2){
			global_density_ub += this->cur_edge_qweight[1 - this->anchor_left].qweight;
			global_density_ub += this->verQWeights[i].qweight;
			if(global_density_ub < this->Minimum_density_needed){
				if(global_density_ub > this->cur_maximum_density_ub){
					this->cur_maximum_density_ub = global_density_ub;
				}
				continue;
			}
		}
		else{break;}
#endif
		
#ifdef SORT_BY_QWEIGHTS
#ifdef GLOBAL_PRUNE
		//cout << "enter global pruning" << endl;
		for(unsigned j = i; m < (this->t - 1) && j < this->n_vertices_local[1 - this->anchor_left]; j++,m++){
			global_density_ub += this->verQWeights[j].qweight;
		}
		if(m == this->t - 1){
			global_density_ub += this->cur_edge_qweight[1 - this->anchor_left].qweight;
			if(global_density_ub < this->Minimum_density_needed){
				if(global_density_ub > this->cur_maximum_density_ub){
					this->cur_maximum_density_ub = global_density_ub;
				}
				break;
			}
		}
		else{
			break;
		}
#endif
#endif
		/**
		 * local pruning
		 *
		 */
		//step 1 : clear L,R,Cand_of_L
		this->L.resize(this->t - 1);
		this->R.resize(this->t - 1);
		this->Cand_of_L.resize(this->t - 1);
		this->aux_vec_R.resize(this->t - 1);
		this->candidate_one_hop_neighbor.resize(this->t - 1);
		this->candidate_verQWeights.resize(this->t - 1);
		this->candidate_two_hop_adj.resize(this->t - 1);
		this->candidate_vertices_rank.resize(this->t - 1);

		for(unsigned i = 0; i < this->L.size(); i++){
			this->L[i].clear();
			for(unsigned j = 0; j < R[i].size(); j++){
				this->R[i][j].clear();
			}
			this->R[i].clear();
			this->aux_vec_R[i].clear();
			this->Cand_of_L[i].clear();
			this->candidate_one_hop_neighbor[i].clear();
			this->candidate_verQWeights[i].clear();
			this->candidate_two_hop_adj[i].clear();
			this->candidate_vertices_rank[i].clear();
		}
		this->L.clear();
		this->R.clear();
		this->aux_vec_R.clear();
		this->Cand_of_L.clear();
		this->candidate_one_hop_neighbor.clear();
		this->candidate_verQWeights.clear();
		this->candidate_two_hop_adj.clear();
		this->candidate_vertices_rank.clear();

		//step 2 : initial
		unsigned root_vertex = this->sorted_vertices[i];
		this->L[0].push_back(root_vertex);
		vector<Neighbor> root_verx_neig = this->local_one_hop_neighbor[1 - this->anchor_left][root_vertex];
		this->R[0].resize(root_verx_neig.size());//num of neighbor
		for(unsigned u = 0; u < root_verx_neig.size(); u++){
			this->R[0][u].push_back(root_verx_neig[u]);
			SmallBclique sbclique(root_verx_neig[u].op_vertex, root_verx_neig[u].related_edge->weight, root_verx_neig[u].related_edge->timestamp);
			this->aux_vec_R[0].push_back(sbclique);
		}
		this->Cand_of_L[0] = this->two_hop_adj_map[root_vertex];
		//step 3 : start seraching
		this->PQbclique(1);

	}

	//step 4: update density upper bound of currently edge
	if(this->insert_or_delete){
		if(!this->cur_maximum_density_ub){
			//cout << "The maximum density upper bound of the current edge is 0" << endl;
			return ;
		}
		this->map1_of_densityUb.insert({this->cur_edge, this->cur_maximum_density_ub});
		this->map2_of_densityUb.insert({this->cur_maximum_density_ub, this->cur_edge});
	}
	else{
		Wtype old_maximum_density_ub = this->map2_of_densityUb.begin()->first;
		if(this->cur_maximum_density_ub < old_maximum_density_ub){
			if(!this->cur_maximum_density_ub){
				auto itr1 = this->map1_of_densityUb.find(this->cur_edge);
				if(itr1 != this->map1_of_densityUb.end()){
					this->map1_of_densityUb.erase(this->cur_edge);
				}
				this->map2_of_densityUb.erase(this->map2_of_densityUb.begin());
			}
			else{
				this->map1_of_densityUb[this->cur_edge] = this->cur_maximum_density_ub;
				this->map2_of_densityUb.erase(this->map2_of_densityUb.begin());
				this->map2_of_densityUb.insert({this->cur_maximum_density_ub, this->cur_edge});
			}
		}
	}
}

void topkbc::PQbclique(unsigned l){
	if(l == this->t - 1){
		unsigned aux_vec_R_size = this->aux_vec_R[l-1].size();
		if(aux_vec_R_size < this->op_t){
			return;
		}
		vector< vector<Neighbor> > temp_R;
		vector<SmallBclique> temp_aux_R;
		SmallBclique op_vert_sbclique;
		vector<Neighbor> op_vert_Neighbors;
		for(unsigned i = 0; i < aux_vec_R_size; i++){
			if(this->aux_vec_R[l-1][i].R_vertex != this->op_vert){
				temp_R.emplace_back(this->R[l-1][i]);
				temp_aux_R.emplace_back(this->aux_vec_R[l-1][i]);
			}
			else{
				op_vert_sbclique = this->aux_vec_R[l-1][i];
				op_vert_Neighbors = this->R[l-1][i];
			}
		}
		unsigned comb_size = this->op_t - 1;
		vector< vector<Neighbor> > temp_R_comb(comb_size);
		vector<SmallBclique> temp_aux_R_comb(comb_size);
		this->find_combinations(temp_R, temp_aux_R, temp_R_comb, temp_aux_R_comb, op_vert_sbclique, op_vert_Neighbors, 0, comb_size, comb_size);

	}
	else{
		if(!this->construct_local_graph(l, this->aux_vec_R[l-1], this->Cand_of_L[l-1])){
			return;
		}
#ifdef LOCAL_PRUNE
		//cout << "enter local pruning" << endl;
		Wtype local_density_ub = 0;
		Wtype ub1 = 0;
		Wtype ub2 = 0;
		
		//step 1 :compute ub1 and ub3 of local density
		vector<SmallBclique> tmp_aux_R;
		tmp_aux_R.assign(this->aux_vec_R[l-1].begin(), aux_vec_R[l-1].end());
		unsigned aux_R_size = tmp_aux_R.size();
		for(unsigned i = 0; i < aux_R_size; i++){
			if(this->cur_edge_one_hop_neighbor[1 - this->anchor_left].find(tmp_aux_R[i].R_vertex) != this->cur_edge_one_hop_neighbor[1 - this->anchor_left].end()){
				tmp_aux_R[i].Lweight += this->cur_edge_one_hop_neighbor[1 - this->anchor_left][tmp_aux_R[i].R_vertex]->weight;
			}
			else{
				tmp_aux_R[i].Lweight += this->cur_edge->weight;
			}
		}
		sort(tmp_aux_R.begin(), tmp_aux_R.end(), cmpByLWeightsDec);
		unsigned cnt = 0;
		bool flag = false;
		for(unsigned i = 0; i < aux_R_size; i++){
			if(cnt < this->op_t - 1){
				ub1 += tmp_aux_R[i].Lweight;
				cnt++;
				if(this->op_vert == tmp_aux_R[i].R_vertex){
					flag = true;
				}
			}
			else{
				if(flag){
					ub1 += tmp_aux_R[i].Lweight;
					break;
				}
				else{
					if(this->op_vert != tmp_aux_R[i].R_vertex){
						continue;
					}
					else{
						ub1 += tmp_aux_R[i].Lweight;
						break;
					}
				}
			}
		}

		//step 2 :compute ub2 of local density
		for(unsigned i = 0; i < this->t-1-l; i++){
			ub2 += this->candidate_verQWeights[l-1][i].qweight;
		}

		//step 3 :compute local density
		local_density_ub = ub1 + ub2;

		if(local_density_ub < this->Minimum_density_needed){//first update maximum density ub
			if(local_density_ub > this->cur_maximum_density_ub){
				this->cur_maximum_density_ub = local_density_ub;
			}
			return;
		}
#endif
		//continue to serach
		unsigned n_cand = this->Cand_of_L[l-1].size();
		for(unsigned i = 0; i < n_cand; i++){
			/**
			 * update L, R, aux_vec_R, C
			 */
			if(i == n_cand + 2 + l - this->t){
				return;
			}
#ifdef SORT_BY_QWEIGHTS_DY
#ifdef LOCAL_PRUNE
			//cout << "enter local pruning" << endl;
			if(i > 0){
				if(i+t-1-l > n_cand){cerr << "error-(i+t-1-l > n_cand)" << endl;getchar();}
				ub2 = 0;
				for(unsigned j = i; j < i+t-1-l; j++){
					ub2 += this->candidate_verQWeights[l-1][j].qweight;
				}
				local_density_ub = ub1 + ub2;
				if(local_density_ub < this->Minimum_density_needed){
					if(local_density_ub > this->cur_maximum_density_ub){
						this->cur_maximum_density_ub = local_density_ub;
					}
					return;
				}
			}
#endif
#endif
			this->L[l] = this->L[l - 1];
			this->L[l].push_back(this->Cand_of_L[l - 1][i]);

			this->R[l] = this->intersection(l, this->R[l-1], this->candidate_one_hop_neighbor[l-1][this->Cand_of_L[l-1][i]], this->aux_vec_R[l-1], 0, 0);
			//cout << "this->aux_vec_R[l].size() = " << this->aux_vec_R[l].size() << endl;
			if(this->aux_vec_R[l].size() < this->op_t){
				continue;
			}
			if(this->candidate_two_hop_adj[l-1].find(this->Cand_of_L[l - 1][i]) != this->candidate_two_hop_adj[l-1].end()){
				this->Cand_of_L[l] = this->candidate_two_hop_adj[l-1][this->Cand_of_L[l - 1][i]];
			}
			else{
				this->Cand_of_L[l].clear();
			}
			//cout << "this->Cand_of_L[l] size = " << this->Cand_of_L[l].size() << " this->L[l].size() = " << this->L[l].size() <<endl;
			if(this->L[l].size() + this->Cand_of_L[l].size() >= this->t - 1){
				//cout << "l + 1 = " << l+1 << " this->t = " << this->t << endl;
				this->PQbclique(l + 1);
			}
		}
	}
}

void topkbc::find_combinations(vector< vector<Neighbor> > &temp_R, vector<SmallBclique> &temp_aux_R, vector< vector<Neighbor> > &temp_R_comb,
		vector<SmallBclique> &temp_aux_R_comb, SmallBclique &op_vert_sbclique,vector<Neighbor> &op_vert_Neighbors,
		unsigned beg_offset, unsigned curr_depth, unsigned total_depth){
    //comb用于临时存储结果。len(comb)==total_depth；beg_offset为左侧游标，初始值取0；
    // total_depth是取出个数；curr_depth用于指示递归深度，初始值取total_depth）
    unsigned N = temp_R.size();
    //cout << "N = " << temp_R.size() << endl;
    if (curr_depth == 0) {
    	Wtype density = this->cur_edge->weight;
     	lint min_ts = this->cur_edge->timestamp;
     	vector<Edge*> edges_of_bclique;

     	edges_of_bclique.emplace_back(this->cur_edge);
    	for(unsigned i = 0; i < total_depth; i++){
    		density += this->cur_edge_one_hop_neighbor[1-this->anchor_left][temp_aux_R_comb[i].R_vertex]->weight;
    		edges_of_bclique.emplace_back(this->cur_edge_one_hop_neighbor[1-this->anchor_left][temp_aux_R_comb[i].R_vertex]);
    		min_ts = min(this->cur_edge_one_hop_neighbor[1-this->anchor_left][temp_aux_R_comb[i].R_vertex]->timestamp, min_ts);

    		//cout << "R_vertex : " << temp_aux_R_comb[i].R_vertex <<endl;
    		//cout << "weight = " << this->cur_edge_one_hop_neighbor[1-this->anchor_left][temp_aux_R_comb[i].R_vertex]->weight << endl;
    		density += temp_aux_R_comb[i].Lweight;
    		//cout << "Lweight : " << temp_aux_R_comb[i].Lweight <<endl;
    		//cout << "temp_aux_R_comb[i].min_ts = " << temp_aux_R_comb[i].min_ts <<endl;
    		if(min_ts > temp_aux_R_comb[i].min_ts){
    			min_ts = temp_aux_R_comb[i].min_ts;
    		}
    		for(unsigned j = 0; j < temp_R_comb[i].size(); j++){
    			edges_of_bclique.emplace_back(temp_R_comb[i][j].related_edge);
    		}
    	}
    	//cout << "op_vert_sbclique Lweight = " << op_vert_sbclique.Lweight << endl;
    	density += op_vert_sbclique.Lweight;

    	if(min_ts > op_vert_sbclique.min_ts){min_ts = op_vert_sbclique.min_ts;}
    	//sencond update maximum density ub
    	if(density < this->Minimum_density_needed){
    		if(density > this->cur_maximum_density_ub){
    			this->cur_maximum_density_ub = density;
    		}
    		return;
    	}
    	else{
    		//if(cur_edge->timestamp == 875726786)
    		//cout << "#find a top-k biclique->density = " << density << " ";
    		if(!this->insert_or_delete){
    			Wtype density_threshold = this->map2_of_densityUb.begin()->first;
    			//cout << "density_threshold = " << density_threshold << endl;
    			if(density > density_threshold){
    				return;
    			}
    		}

    		for(unsigned i = 0; i < op_vert_Neighbors.size(); i++){
    			edges_of_bclique.emplace_back(op_vert_Neighbors[i].related_edge);
    			min_ts = min(op_vert_Neighbors[i].related_edge->timestamp, min_ts);
    		}
    		//cout << "min_ts = " << min_ts << endl;
    		pqbclique pqb(this->cur_edge, min_ts, edges_of_bclique);
    		this->map_of_Topk_res.insert(make_pair(density, pqb));
    		if(this->map_of_Topk_res.size() == this->k){
    			this->Minimum_density_needed = (*this->map_of_Topk_res.rbegin()).first;
			}
    		if(this->map_of_Topk_res.size() > this->k){
    			Wtype replaced_val_of_curEdge = this->replace_Topk_and_update_densityUbMap();
    			//last update maximum density ub
    			if(replaced_val_of_curEdge > this->cur_maximum_density_ub){
					this->cur_maximum_density_ub = replaced_val_of_curEdge;
				}
    			this->Minimum_density_needed = (*this->map_of_Topk_res.rbegin()).first;
    		}
    	}
        return;
    }
    for (unsigned i = beg_offset; i < N; i++){
    	temp_R_comb[total_depth-curr_depth] = temp_R[i];
    	temp_aux_R_comb[total_depth-curr_depth] = temp_aux_R[i];
        //comb[total_depth-curr_depth] = seed_vertices[i];
        find_combinations(temp_R, temp_aux_R, temp_R_comb, temp_aux_R_comb, op_vert_sbclique, op_vert_Neighbors, i + 1, curr_depth - 1, total_depth);
    }
}

Wtype topkbc::replace_Topk_and_update_densityUbMap(){
	Wtype kth_value = 0;
	Wtype expire_value = 0;
	Wtype replaced_val_of_curEdge = 0;
	int offset = this->map_of_Topk_res.size() - this->k;
	bool is_replace = 0;

	//step 1:  replace Topk
	for(auto itr1 = this->map_of_Topk_res.rbegin(); itr1 != this->map_of_Topk_res.rend();){
		if(offset == 1){
			expire_value = itr1->first;
		}
		itr1 ++;
		offset--;
		if(offset == 0){
			kth_value = itr1->first;
			if(kth_value == expire_value){
				return 0;
			}
			else{
				is_replace = 1;
				//cout << "expire value = " <<  expire_value << endl;
			}
		}
	}

	//step 2: update map of density upper bound
	//cout << "1-map1_of_densityUb size = " <<this->map1_of_densityUb.size() << endl;
	//cout << "1-map2_of_densityUb size = " <<this->map2_of_densityUb.size() << endl;
	if(is_replace){
		for(;;){
			auto itr1 = this->map_of_Topk_res.rbegin();
			if(itr1->first <= expire_value){
				if(itr1->second.related_edge != this->cur_edge){
					if(this->map1_of_densityUb.find(itr1->second.related_edge) != this->map1_of_densityUb.end()){
						Wtype density_ub = this->map1_of_densityUb[itr1->second.related_edge];
						this->map1_of_densityUb[itr1->second.related_edge] = itr1->first;
						auto ret = this->map2_of_densityUb.equal_range(density_ub);
						auto it = ret.first;
						while(it != ret.second){
							if(it->second == itr1->second.related_edge){
								this->map2_of_densityUb.erase(it);
								this->map2_of_densityUb.insert({itr1->first, itr1->second.related_edge});
								break;
							}
							else{
								it++;
//								cout << "error!!! not found in map2_of_densityUb" <<endl;
//								getchar();
							}
						}
					}
					else{
						this->map1_of_densityUb.insert({itr1->second.related_edge, itr1->first});
						this->map2_of_densityUb.insert({itr1->first,itr1->second.related_edge});
						//cout << "not find in map1_of_densityUb" <<endl;
					}
				}
				else{
					replaced_val_of_curEdge = expire_value;
				}
				this->map_of_Topk_res.erase((++itr1).base());
				continue;
			}
			else{
				return replaced_val_of_curEdge;
			}
		}
	}

	//cout << "OUT replace_Topk_and_update_densityUbMap()" << endl;
	return replaced_val_of_curEdge;
}

void topkbc::remove_expired_result(){
	// step 1: remove old Top-k (p,q)-bcliques
	for(auto map_itr = this->map_of_Topk_res.begin(); map_itr != this->map_of_Topk_res.end();){
		//cout << "map_itr->second.min_timestamp = " << map_itr->second.min_timestamp << endl;
		if(map_itr->second.min_timestamp < this->win_start_time){
			this->map_of_Topk_res.erase(map_itr ++);
		}
		else{
			map_itr ++;
		}
	}
	
	// step 2: remove old edge and map of densityub
	for(int i = 0; i < (int)this->snapshot.size();){
		if(this->snapshot[i]->timestamp < this->win_start_time){//如果时间戳与过期边的相同，则删除
			auto map_itr = this->map1_of_densityUb.find(this->snapshot[i]);
			if(map_itr != this->map1_of_densityUb.end()){
				auto ret = this->map2_of_densityUb.equal_range(map_itr->second);
				auto it = ret.first;
				while(it != ret.second){
					if(it->second == map_itr->first){
						this->map2_of_densityUb.erase(it);
						break;
					}
					else{
						it++;
					}
				}
				this->map1_of_densityUb.erase(map_itr);
			}
			this->construct_global_index(this->snapshot[i], 0);
			this->updated_edges_delete ++;
			//this->num_edges--; //compute num of edges

			delete this->snapshot[i];
			this->snapshot[i]=NULL;
			this->snapshot.erase(this->snapshot.begin() + i);
		}
		else{
			break;
		}
	}
	
/**
	this->tmp_ts = this->snapshot[this->stride - 1]->number;

	// step 1: remove old bclique
	for(auto map_itr = this->map_of_Topk_res.begin(); map_itr != this->map_of_Topk_res.end();){
		//cout << "map_itr->second.min_timestamp = " << map_itr->second.min_timestamp << endl;
		if(map_itr->second.min_timestamp <= this->tmp_ts){
			this->map_of_Topk_res.erase(map_itr ++);
		}
		else{
			map_itr ++;
		}
	}

	// step 2: remove old edge and map of densityub
	for(int i = 0; i < (int)this->snapshot.size();){
		if(this->snapshot.size() == 0){
			cerr << "error" << endl;
			system("pause");
		}
		if(this->snapshot[i]->number <= this->tmp_ts){//如果时间戳与过期边的相同，则删除
			auto map_itr = this->map1_of_densityUb.find(this->snapshot[i]);
			if(map_itr != this->map1_of_densityUb.end()){
				auto ret = this->map2_of_densityUb.equal_range(map_itr->second);
				auto it = ret.first;
				while(it != ret.second){
					if(it->second == map_itr->first){
						this->map2_of_densityUb.erase(it);
						break;
					}
					else{
						it++;
					}
				}
				this->map1_of_densityUb.erase(map_itr);
			}
			this->construct_global_index(this->snapshot[i], 0);
			//this->num_edges--; //compute num of edges

			delete this->snapshot[i];
			this->snapshot[i]=NULL;
			this->snapshot.erase(this->snapshot.begin() + i);
		}
		else{
			break;
		}
	}
	//cout << "OUT remove_expired_result()" << endl;
	*/
}

void topkbc::search_alternate_biclique(){
	this->insert_or_delete = false;
	if(this->map2_of_densityUb.empty() || this->map2_of_densityUb.begin()->first == 0){
		//cout << "do not need to search alternate biclique" << endl;
		return;
	}
	Edge* serach_in_edge = nullptr;
	this->next_maximum_density_ub = 0;
	for(;;){
		auto map_itr = this->map_of_Topk_res.rbegin();
		if((this->map_of_Topk_res.size() >= this->k && map_itr->first > this->next_maximum_density_ub)|| this->map2_of_densityUb.empty()){
			return;
		}
		auto temp_itr = this->map2_of_densityUb.begin();
		serach_in_edge = temp_itr->second;
		temp_itr ++;
		
		if(temp_itr == this->map2_of_densityUb.end()){
			//cout << "map2_of_densityUb have no next_maximum_density_ub" << endl;
			this->next_maximum_density_ub = 0;
		}
		else{
			this->next_maximum_density_ub = temp_itr->first;
		}
		this->cur_edge = serach_in_edge;
		
		if(!this->construct_local_graph(serach_in_edge, 0)){
			this->map2_of_densityUb.erase(this->map2_of_densityUb.begin());
			auto temp_itr1 = this->map1_of_densityUb.find(serach_in_edge);
			if(temp_itr1 != this->map1_of_densityUb.end()){
				this->map1_of_densityUb.erase(temp_itr1);
			}
			continue;
		}
		else{
			this->serach_Topk_bcliques();
		}

	}
	//cout << "OUT search_alternate_biclique()" << endl;
}

void topkbc::print_local_one_hop_neighbor(){
	cout << "#Local graph left vertices num = " << this->n_vertices_local[0] << ", list as follow :" << endl;
	for(auto itr1 : this->local_one_hop_neighbor[0]){
		cout << itr1.first << " | ( ";
		for(auto itr2 : itr1.second){
			cout << itr2.op_vertex << " ";
		}
		cout << ")" << endl;
	}

	cout << "#Local graph right vertices num = " << this->n_vertices_local[1] << endl;
	for(auto itr1 : this->local_one_hop_neighbor[1]){
		cout << itr1.first << " | ( ";
		for(auto itr2 : itr1.second){
			cout << itr2.op_vertex << " ";
		}
		cout << ")" << endl;
	}
	cout << endl;
}

void topkbc::print_sorted_vertices(){
	cout << "#sorted vertices = ( ";
	for(unsigned i = 0; i < this->verQWeights.size(); i++){
		cout << this->verQWeights[i].vertex << "|" << this->verQWeights[i].qweight << " ";
	}
	cout << ")" << endl;
	cout << endl;
}

void topkbc::print_local_two_hop_neighbor(){
	cout << "#two hop neighbor listed as follow : " << endl;
	for(auto itr1 : this->two_hop_adj_map){
		cout << "( " << itr1.first << " | ";
		for(unsigned i = 0; i < itr1.second.size(); i++){
			cout << itr1.second[i] << " ";
		}
		cout << ")" <<endl;
	}
}

void topkbc::print_cur_edge_one_hop_neighbor(){
	cout <<"#currently edge one hop neighbor list as follow: " <<endl;
	cout << "u = " << this->cur_edge->u << " v = " << this->cur_edge->v <<endl;
	for(auto itr : this->cur_edge_one_hop_neighbor[0]){
		cout << itr.first << "-->weight: " << itr.second->weight <<endl;
	}
	cout << "----------------------------" << endl;
	for(auto itr : this->cur_edge_one_hop_neighbor[1]){
		cout << itr.first << "-->weight: " << itr.second->weight <<endl;
	}
	cout << "----------------------------" << endl;
}

void topkbc::print_Topk_weight_bicliques(){
	cout << "Total Top-k bicliques# results : " << this->map_of_Topk_res.size() << endl;
	cout << "The current updated number of edges is : " << this->num_edges << endl;
	cout << "Initial runing time: using " << this->time_initial << " seconds." << endl;
	//cout << "---------------------" << endl;
#ifdef TOPK_RESULT_OUTPUT_FILE
	//cout << "Total Top-k bicliques# results : " << this->map_of_Topk_res.size() << endl;
	unsigned u = 0;
	unsigned v = 0;
	vector<unsigned> l_vertices;
	vector<unsigned> r_vertices;
	string dataname = data_file_path.substr(0, data_file_path.rfind("."));
	cout << "dataname : " << dataname << endl;
	string out_file_path = this->make_outfile_name(dataname, p, q, this->update_count, this->map_of_Topk_res.size(), this->winsz);
	cout << out_file_path << endl;
	ofstream fout;
	fout.open(out_file_path);
	for(auto it1 : this->map_of_Topk_res){
		l_vertices.clear();
		r_vertices.clear();
		fout << "density: " << it1.first << " || " << "vertices: ";
		for(unsigned i = 0; i < it1.second.edges_of_bclique.size(); i++){
			u = it1.second.edges_of_bclique[i]->u;
			v = it1.second.edges_of_bclique[i]->v;
			//cout << "(" << u << " " << v << ")" << endl;
			if(find(l_vertices.begin(), l_vertices.end(), u) == l_vertices.end()){l_vertices.push_back(u);}
			if(find(r_vertices.begin(), r_vertices.end(), v) == r_vertices.end()){r_vertices.push_back(v);}
		}
		sort(l_vertices.begin(), l_vertices.end());
		sort(r_vertices.begin(), r_vertices.end());
		for(unsigned i = 0; i < l_vertices.size(); i++){
			fout << l_vertices[i] << " ";
		}
		fout << "| ";
		for(unsigned i = 0; i < r_vertices.size(); i++){
			fout << r_vertices[i] << " ";
		}
		fout << endl;
	}
	fout.close();
#endif
}

void topkbc::print_map_of_edgeRank(){
	
}

void topkbc::print_edges(){
//	string dataname = data_file_path.substr(0, data_file_path.rfind("."));
//	cout << "dataname : " << dataname << endl;
//	string out_file_path = this->make_outfile_name(dataname, p, q, this->update_count, this->map_of_edges.size(), this->num_edges);
//	cout << out_file_path << endl;
//	ofstream fout;
//	fout.open(out_file_path);
//	for(auto it : this->map_of_edges){
//		fout << it.first.first << " " << it.first.second << endl;
//	}
//	fout.close();
}



string topkbc::make_outfile_name(string dataname, unsigned pv, unsigned qv, lint update_num, unsigned clique_num, unsigned edge_num){
    string ans = "";
    
    std::size_t found = dataname.find_last_of("/");
#ifdef SORT_BY_QWEIGHTS_DY
    ans += dataname.substr(0, found) + "/result/our-final/" + dataname.substr(found + 1) + "_" + to_string(pv) + "_" + to_string(qv) + "_" + to_string(update_num) + "_" + to_string(clique_num) + "_" + to_string(edge_num) + ".txt";
    return ans;
#endif
#ifdef SORT_BY_QWEIGHTS
    ans += dataname.substr(0, found) + "/result_bckeep_qweight/" + dataname.substr(found + 1) + "_" + to_string(pv) + "_" + to_string(qv) + "_" + to_string(update_num) + "_" + to_string(clique_num) + "_" + to_string(edge_num) + ".txt";
#endif
#ifdef SORT_BY_DEGREE
    ans += dataname.substr(0, found) + "/result_bckeep_degree/" + dataname.substr(found + 1) + "_" + to_string(pv) + "_" + to_string(qv) + "_" + to_string(update_num) + "_" + to_string(clique_num) + "_" + to_string(edge_num) + ".txt";
#endif
#ifdef SORT_BY_CORE
    ans += dataname.substr(0, found) + "/result_bckeep_core/" + dataname.substr(found + 1) + "_" + to_string(pv) + "_" + to_string(qv) + "_" + to_string(update_num) + "_" + to_string(clique_num) + "_" + to_string(edge_num) + ".txt";
#endif
    return ans;
}

vector<unsigned> topkbc::intersection(vector<unsigned> &vec1, vector<unsigned> &vec2, int offset1, int offset2){
    vector<unsigned> ans;
    if(vec1.empty() || vec2.empty() || offset1 >= vec1.size() || offset2 >= vec2.size()){
        return ans;
    }

    // we assume vec1 and vec2 are sorted by ascending order of the elements
    unsigned i, j;
    i = offset1;
    j = offset2;
    while(i < vec1.size() && j < vec2.size()){
        if(vec1[i] == vec2[j]){
            ans.emplace_back(vec1[i]);
            i++;
            j++;
        }
        else if(vec1[i] < vec2[j]){
            i++;
        }
        else{
            j++;
        }
    }
    return ans;
}

vector<Neighbor> topkbc::intersection(bool insert_or_del, vector<Neighbor> &vec1, vector<Neighbor> &vec2, int offset1, int offset2){
    vector<Neighbor> ans;
    if(vec1.empty() || vec2.empty() || offset1 >= vec1.size() || offset2 >= vec2.size()){
        return ans;
    }

    // we assume vec1 and vec2 are sorted by ascending order of the elements
    unsigned i, j;
    i = offset1;
    j = offset2;
    if(insert_or_del){
		while(i < vec1.size() && j < vec2.size()){
			if(vec1[i] == vec2[j]){
				ans.emplace_back(vec1[i]);
				i++;
				j++;
			}
			else if(vec1[i] < vec2[j]){
				i++;
			}
			else{
				j++;
			}
		}
    }
    else{
    	while(i < vec1.size() && j < vec2.size()){
    		if(vec1[i] == vec2[j] && vec1[i].related_edge->number <= this->cur_edge->number && vec2[j].related_edge->number <= this->cur_edge->number){
    			ans.emplace_back(vec1[i]);
    			//cout << vec1[i].op_vertex << " ";
    			i++;
    			j++;
    		}
    		else if(vec1[i] < vec2[j]){
				i++;
			}
			else{
				j++;
			}
    	}
    	//cout << endl;
    }
    return ans;
}

vector<vector<Neighbor>> topkbc::intersection(unsigned l, vector<vector<Neighbor>> &vec1, vector<Neighbor> &vec2, vector<SmallBclique> &vec3, int offset1, int offset2){
	vector<vector<Neighbor>> ans;
	vector<Neighbor> same_neig;
	vector<SmallBclique> same_sbclique;
	unsigned R_vertex;
	lint min_ts;
	Wtype Lweight;
    if(vec1.empty() || vec2.empty() || offset1 >= vec1.size() || offset2 >= vec2.size()){
        return ans;
    }

    // we assume vec1 and vec2 are sorted by ascending order of the elements
    unsigned i, j;
    i = offset1;
    j = offset2;
    while(i < vec1.size() && j < vec2.size()){
        if(vec1[i][0] == vec2[j]){
        	same_neig = vec1[i];
        	same_neig.push_back(vec2[j]);
            ans.emplace_back(same_neig);
            R_vertex = vec2[j].op_vertex;
            Lweight = vec3[i].Lweight + vec2[j].related_edge->weight;
            min_ts = min(vec2[j].related_edge->timestamp, vec3[i].min_ts);
            SmallBclique sbc(R_vertex, Lweight, min_ts);
            same_sbclique.push_back(sbc);
            i++;
            j++;
        }
        else if(vec1[i][0] < vec2[j]){
            i++;
        }
        else{
            j++;
        }
    }
    this->aux_vec_R[l] = same_sbclique;
    return ans;
}

vector<Neighbor> topkbc::intersection(vector<SmallBclique> &vec1, vector<Neighbor> &vec2, map<unsigned, vector<unsigned>> &R_one_hop_neig, unsigned vertex, int offset1, int offset2){
	vector<Neighbor> ans;
	if(vec1.empty() || vec2.empty() || offset1 >= vec1.size() || offset2 >= vec2.size()){
		return ans;
	}

	// we assume vec1 and vec2 are sorted by ascending order of the elements
	unsigned i, j;
	i = offset1;
	j = offset2;
	while(i < vec1.size() && j < vec2.size()){
		if(vec1[i].R_vertex == vec2[j].op_vertex){
			R_one_hop_neig[vec1[i].R_vertex].push_back(vertex);
			ans.emplace_back(vec2[j]);
			i++;
			j++;
		}
		else if(vec1[i].R_vertex < vec2[j].op_vertex){
			i++;
		}
		else{
			j++;
		}
	}
	return ans;
}

vector<Neighbor> topkbc::merge(vector<Neighbor> &vec1, vector<Neighbor> &vec2){
    vector<Neighbor> ans;

    // we assume vec1 and vec2 are sorted by ascending order of the elements
    unsigned i, j;
    i = 0;
    j = 0;
    while(i < vec1.size() && j < vec2.size()){
        if(vec1[i] == vec2[j]){
            ans.emplace_back(vec1[i]);
            i++;
            j++;
        }
        else if(vec1[i] < vec2[j]){
            ans.emplace_back(vec1[i]);
            i++;
        }
        else{
            ans.emplace_back(vec2[j]);
            j++;
        }
    }
    return ans;
}

unsigned topkbc::intersection_count(vector<unsigned> &vec1, vector<unsigned> &vec2){
    unsigned ans = 0;
    if(vec1.empty() || vec2.empty()){
        return 0;
    }

    // we assume vec1 and vec2 are sorted by ascending order of the elements
    unsigned i, j;
    i = j = 0;
    while(i < vec1.size() && j < vec2.size()){
        if(vec1[i] == vec2[j]){
            ans++;
            i++;
            j++;
        }
        else if(vec1[i] > vec2[j]){
            j++;
        }
        else{
            i++;
        }
    }

    return ans;
}

long topkbc::getMemoryUse(){
    int who = RUSAGE_SELF;
    struct rusage usage;
    getrusage(who, &usage);
    return usage.ru_maxrss;
}