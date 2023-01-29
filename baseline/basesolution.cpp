#include "basesolution.h"


basesolution::basesolution(int _p, int _q, int _winsz, string _data_path, long long _update_times, int _k)
{
		this->p = _p;
		this->q = _q;
		this->winsz = _winsz;
		this->data_path = _data_path;
		this->update_times =_update_times;
		this->update_cnt = 0;
		this->k = _k;
		this->ts_is_same = 0;
		this->tmp_ts = 0;

		for(int i = 0; i < this->p; i ++){
			this->r[i] = new vector<int>;
			this->c[i] = new vector<int>;
		}
		for(int i = this->p; i < MAX_P; i ++){
			this->r[i] = NULL;
			this->c[i] = NULL;
		}
}

basesolution::~basesolution(){
	for(int i = 0; i < (int)this->snapshot.size();i ++){
		delete this->snapshot[i];
		this->snapshot[i] = NULL;
	}
	this->snapshot.clear();
	for(int i = 0; i < MAX_P; i ++){
		delete this->r[i];
		this->r[i] = nullptr;
		delete this->c[i];
		this->c[i] = nullptr;
	}

}

void basesolution::run(){
	gstream gs(this->data_path);
	Edge* newedge = nullptr;
	while(gs.hasnext()){
		newedge = gs.next();
		if(this->update_cnt >= this->update_times && newedge->timestamp != this->snapshot[this->snapshot_size()-1]->timestamp){
			cout << "this->Resultmap.size() is " << this->Resultmap.size() << endl;
			cout << "update times already get Maxvalue = " << this->update_cnt <<endl;
			break;
		}
		if(this->snapshot_size() == 0){
#ifdef DEBUG
			cout<<"init snapshot"<<endl;
#endif
			this->snapshot.push_back(newedge);
			this->update_cnt ++;
			this->update(newedge);
			continue;
		}
#ifdef DEBUG
		cout << "snapshot size = " <<this->snapshot_size() << endl;
		cout << "comming newedge->u = " << newedge->u << endl;
		cout << "comming newedge->v = " << newedge->v << endl;
		cout << "comming newedge->weight = " << newedge->weight << endl;
		//getchar();
#endif
		if(newedge->timestamp == this->snapshot[this->snapshot_size()-1]->timestamp)
		{//同一时刻可能有多条边出现
			this->ts_is_same = 1;
			this->snapshot.push_back(newedge);
			this->update(newedge);
		}
		else{
			this->ts_is_same = 0;
			this->snapshot.push_back(newedge);
			this->update_cnt ++;
			this->update(newedge);
#ifdef DEBUG

			cout << "snapshot size = " <<this->snapshot_size() << endl;
			//getchar();
#endif
		}

	}
}

void basesolution::update(Edge* _newedge){
	if(this->update_cnt <= this->winsz){//更新次数小于等于窗口大小，不考虑过期边，只关注插入
		this->insert(_newedge);
	}
	else{//更新次数大于窗口大小,考虑过期边，同时关注插入和删除
		if(this->ts_is_same){
			this->insert(_newedge);
		}
		else{
			this->remove();
			this->insert(_newedge);
		}
	}
}

void basesolution::insert(Edge* _newedge){
#ifdef DEBUG
	if(this->snapshot_size() == 0){
		cerr << "snapshot size == 0,error,in basesolution::insert" << endl;
		system("pause");
	}
#endif
	//this->localsubg.clear();
	this->neighbor_of_u.clear();
	this->neighbor_of_v.clear();
	this->one_hop_tb.clear();
	//generate neighbor of u and v
	/*
	 * neighbor of v is 2-hop neighbor of u
	 * neighbor of u is superset of R
	 */
	for(int i = 0; i < this->snapshot_size(); i ++){
		if(_newedge->u == this->snapshot[i]->u){
			this->one_hop_tb[_newedge->u].insert(this->snapshot[i]);
			this->neighbor_of_u.insert(this->snapshot[i]->v);
		}
		if(_newedge->v == this->snapshot[i]->v && _newedge->u != this->snapshot[i]->u){
			this->neighbor_of_v.insert(this->snapshot[i]->u);
		}
	}
#ifdef DEBUG
	cout << "this->neighbor_of_u.size() is " << this->neighbor_of_u.size() << endl;
	cout << "this->neighbor_of_v.size() is " << this->neighbor_of_v.size() << endl;
//	if(this->localsubg.size() == 0){
//		cerr << "localsubg size == 0, error" << endl;
//		system("pause");
//	}
//	cout << "this->localsubg.size() is " << this->localsubg.size() <<endl;
#endif
	this->collect_2hop_neighbors(_newedge);
	this->init_partial_biclique(_newedge);
	this->recursionTree_construction(1, this->L);
#ifdef DEBUG
	cout << "this->commom_neighbor_cnt.size() is " << this->commom_neighbor_cnt.size() << endl;
	cout << "this->update_cnt is " << this->update_cnt << endl;
	cout << "this->Edgemap.size() is " << this->Edgemap.size() << endl;
	cout << "edgemap weight as follow :" << endl;
	for(auto itr : this->Edgemap){
		cout << itr.first << endl;
	}
	cout << "this->Resultmap.size() is " << this->Resultmap.size() << endl;
	cout << "The top-k result weight as follow :" << endl;
	for(auto itr : this->Resultmap){
		cout << itr.first << endl;
	}
#endif
}

void basesolution::remove(){
#ifdef DEBUG
	if(this->snapshot_size() == 0){
		cerr << "error" << endl;
		system("pause");
	}
	cout << "start remove()" <<endl;
#endif

	this->tmp_ts = this->snapshot[0]->timestamp;
#ifdef DEBUG
	cout << "this->tmp_ts = " << this->tmp_ts <<endl;
#endif
	this->rm_old_result(this->tmp_ts);//delete expired bicliques
	for(int i = 0; i < this->snapshot_size();){
		if(this->snapshot_size() == 0){
			cerr << "error" << endl;
			system("pause");
		}
		if(this->snapshot[i]->timestamp == this->tmp_ts){//如果时间戳与过期边的相同，则删除
#ifdef DEBUG
			cout << "removeing this->tmp_ts = " << this->tmp_ts << endl;
			cout << "removeing this->snapshot[i]->timestamp = " << this->snapshot[i]->timestamp << endl;
#endif
			delete this->snapshot[i];
			this->snapshot[i]=NULL;
			this->snapshot.erase(this->snapshot.begin() + i);
			//cout << "removing finished success" <<endl;
		}
		else{
			break;
		}
	}
#ifdef DEBUG
	for(auto it : this->snapshot){
		cout << it->u << "--->" << it->v <<endl;
	}
#endif

}

void basesolution::collect_2hop_neighbors(Edge* _newedge){
	this->commom_neighbor_cnt.clear();
	for(auto itr1 : this->neighbor_of_v){
		for(int i = 0; i < (int)this->snapshot_size(); i ++){
			if(itr1 == this->snapshot[i]->u){
				//2-hop neigbor of u neighbor is in neigbor of u
				if(this->neighbor_of_u.find(this->snapshot[i]->v) != this->neighbor_of_u.end()){
					//compute neigbor of all 2-hop neigbor of u
					this->one_hop_tb[itr1].insert(this->snapshot[i]);
					if(this->commom_neighbor_cnt.find(itr1) != this->commom_neighbor_cnt.end()){
						this->commom_neighbor_cnt[itr1] ++;
					}
					else{
						this->commom_neighbor_cnt[itr1] = 1;
					}
				}
			}
		}
	}
	//compute localsubg all vertex one hop neighbor
//	for(int i = 0; i < (int)_localsubg.size(); i ++){
//		this->one_hop_tb[_localsubg[i]->u].insert(_localsubg[i]);
//	}

	//compute two hop neighbor
	//set<Edge*>::iterator itr;
//		for(auto itr1 : this->one_hop_tb[_newedge->u]){
//		for(auto itr2 : this->one_hop_tb[itr1->v]){
//			if(itr2->u != _newedge->u){
//				if(this->commom_neighbor_cnt.find(itr2->u) != this->commom_neighbor_cnt.end()){
//					this->commom_neighbor_cnt[itr2->u] ++;
//				}
//				else{
//					this->commom_neighbor_cnt[itr2->u] = 1;
//				}
//			}
//		}
//	}
	//compute strength >= q or p 2-hop neighbor
	for(auto itr3 = this->commom_neighbor_cnt.begin(); itr3 != this->commom_neighbor_cnt.end();){
		if(itr3->second < this->q){
			this->commom_neighbor_cnt.erase(itr3++);
		}
		else{
			itr3++;
		}
	}
#ifdef DEBUG
	cout << "this is commom_neighbor_cnt" <<endl;
	for(auto it : this->commom_neighbor_cnt){
		cout << it.first<< "--->" << it.second << endl;
	}
	cout << "commom_neighbor_cnt finish traval" <<endl;

#endif
}

void basesolution::init_partial_biclique(Edge* _newedge){
	cout << "In init_partial_biclique " <<endl;
	this->L.clear();
	this->r[0]->clear();
	this->c[0]->clear();
	this->L.push_back(_newedge->u);
	cout << "init R = ( ";
	for(auto it : this->neighbor_of_u){
#ifdef DEBUG
		cout << it << " ";
#endif
		this->r[0]->push_back(it);
	}
	cout << ")" <<endl;
	cout << "init C = (";
	for(auto it : this->commom_neighbor_cnt){
#ifdef DEBUG
		cout << it.first << " ";
#endif
		this->c[0]->push_back(it.first);
	}
	cout << ")" <<endl;
}

void basesolution::neighbor_intersection(int l, int &_c, vector<int>* _R){
	cout << "In neighbor_intersection" <<endl;
	cout << "_c = " << _c <<endl;
	vector<int> tmp_vec;
	_R->clear();
	if(this->one_hop_tb.find(_c) != this->one_hop_tb.end()){
		for(auto itr : this->one_hop_tb[_c]){
			cout << "itr->v =" << itr->v <<endl;
			tmp_vec.push_back(itr->v);
		}
		sort(tmp_vec.begin(), tmp_vec.end());
		sort(this->r[l-1]->begin(), this->r[l-1]->end());
		set_intersection(tmp_vec.begin(), tmp_vec.end(), this->r[l-1]->begin(), this->r[l-1]->end(), back_inserter(*(_R)));
	}
	else{
		cerr<< "void basesolution::neighbor_intersection error" <<endl;
		system("pause");
	}
	cout << "Out neighbor_intersection" <<endl;

}

void basesolution::update_Candidate(int l, int i, vector<int>* _R, vector<int>* _C){
	cout << "In update_Candidate " <<endl;
	vector<int> tmp_vec;
	vector<int> res_vec;
	this->c[l-1+1]->clear();
	cout << "i = " << i <<endl;
	cout << "(int)_C.size() = " << _C->size() << endl;
	if(i == (int)_C->size()){//no candidate can choose
#ifdef DEBUG
		cout << "no candidate can choose" << endl;
#endif
		return;
	}
	else{
		for(; i < (int)_C->size(); i ++){
			for(auto itr : this->one_hop_tb[_C->at(i)]){
				//cout << itr->v << " ";
				tmp_vec.push_back(itr->v);
			}
			cout << "(int)tmp_vec.size() = " << (int)tmp_vec.size() <<endl;
			cout << "_R.size() = " << _R->size() << endl;
			sort(tmp_vec.begin(), tmp_vec.end());
			sort(_R->begin(), _R->end());
			set_intersection(tmp_vec.begin(), tmp_vec.end(), _R->begin(), _R->end(), back_inserter(res_vec));
			//cout << "after set_intersection (int)res_vec.size() = " << (int)res_vec.size() <<endl;
			if((int)res_vec.size() < this->q){
				//cout << "(int)res_vec.size() = " << (int)res_vec.size() <<endl;
				//this->tmp_C.erase(this->tmp_C.begin() + i);
				cout << "_C[i] = " << _C->at(i) << ", res_vec() < q" <<endl;
			}
			else{
				this->c[l-1+1]->push_back(_C->at(i));
				//cout << "tmp_C size = " << this->tmp_C.size() <<endl;
			}
			tmp_vec.clear();
			res_vec.clear();
		}
	}
	cout << "Out update_Candidate " <<endl;
}

void basesolution::recursionTree_construction(int l, vector<int> &_L){
	cout << "recursionTree_construction begin " <<endl;
	cout << " _L.size() = " << _L.size() << endl;
	cout << " r[l-1].size() = " << this->r[l-1]->size() << endl;
	cout << " c[l-1].size() = " << this->c[l-1]->size() << endl;
	if((int)this->r[l-1]->size() < this->q){
		return;
	}
	if(l == this->p){
#ifdef DEBUG
		if(l != (int)_L.size()){
			cerr << "l != _L.size()" <<endl;
			system("pause");
		}
#endif
		this->add_new_result(_L, this->r[l-1]);
		return;
	}
	else{
		for(int i = 0; i < (int)this->c[l-1]->size();){
			if((int)_L.size() + (int)this->c[l-1]->size() < this->p){//
				return;
			}
			this->neighbor_intersection(l, this->c[l-1]->at(i), this->r[l-1+1]);//choose one vertex from _C,e.g._C[i],update tmp_R

			//this->tmp_C.erase(tmp_C.begin()+j);//can not choose previous vertex

			this->update_Candidate(l, i + 1, this->r[l-1+1], this->c[l-1]);//update tmp_C
			if((int)_L.size() + 1 + (int)this->c[l-1+1]->size() < this->p || (int)this->r[l-1+1]->size() < this->q){
				cout << "(int)this->c[l-1+1]->size() = " << (int)this->c[l-1+1]->size() <<endl;
				cout << "(int)this->r[l-1+1]->size() = " << (int)this->r[l-1+1]->size() <<endl;
				this->c[l-1]->erase(this->c[l-1]->begin() + i);//can not choose previous vertex
				continue;
			}
			else{
				_L.push_back(this->c[l-1]->at(i));
				//_R = this->tmp_R;
				//_C = this->tmp_C;
#ifdef DEBUG
				cout << "the next recursion : " <<endl;
				cout << "(int)_L.size() = " << (int)_L.size() <<endl;
				cout << "(int)_R.size() = " << (int)this->r[l-1+1]->size() <<endl;
				cout << "(int)_C.size() = " << (int)this->c[l-1+1]->size() <<endl;
				cout << "l + 1 = " << l + 1 <<endl;
#endif
				recursionTree_construction(l + 1, _L);
				//backtrack
				_L.pop_back();
				i ++;

			}
		}

	}
}

void basesolution::find_combinations(vector< vector<int> > &combs, int cur_picked, int min_subscript, vector<int> &_R, vector<int> &comb){
	if(cur_picked == this->q - 1){
		combs.push_back(comb);
		return;
	}
	for(int i = min_subscript; i <= (int)_R.size() - this->q + 1 + cur_picked; i ++){
		comb[cur_picked] = _R[i];
		find_combinations(combs, cur_picked + 1, i + 1, _R, comb);
	}
}

void basesolution::rm_old_result(long long _timestamp){
#ifdef DEBUG
	cout << "In rm_old_result"  <<endl;
#endif
	if(this->Resultmap.size() == 0 ){
		return ;
	}
	else{
		bool flag = 0;
		for(auto itr1 = this->Resultmap.begin(); itr1 != this->Resultmap.end();){
			cout << "(*itr1).second size is " << (*itr1).second.size() <<endl;
			for(auto itr2 : (*itr1).second){
				cout << "u = " <<itr2->u << ", v = " << itr2->v <<endl;
				cout << "itr2->timestamp = " << itr2->timestamp <<endl;
				if(_timestamp == itr2->timestamp){
					cout << "itr2->timestamp == _timestamp == " << itr2->timestamp << endl;
					this->Resultmap.erase(itr1 ++);
					flag = 1;
					break;
				}
			}
			if(itr1 != this->Resultmap.end()){
				if(!flag){
					itr1 ++;
				}
				else{
					flag = 0;
				}
			}
			else{
				break;
			}
		}
	}
	if(this->Edgemap.size() == 0){
		return ;
	}
	else{
		bool flag = 0;
		for(auto itr1 = this->Edgemap.begin(); itr1 != this->Edgemap.end();){
			for(auto itr2 : (*itr1).second){
				if(_timestamp == itr2->timestamp){
					this->Edgemap.erase(itr1 ++);
					flag = 1;
					break;
				}
			}
			if(itr1 != this->Resultmap.end()){
				if(!flag){
					itr1 ++;
				}
				else{
					flag = 0;
				}
			}
			else{
				break;
			}
		}
		if(this->Edgemap.size() != 0){
			cout << "this->Edgemap.size() > 0" << endl;
			if((int)this->Resultmap.size() < this->k){
				auto itr3 = this->Edgemap.begin();
				while(itr3 != this->Edgemap.end()){
					this->Resultmap.insert(multimap<int,set<Edge*>>::value_type((*itr3).first, (*itr3).second));
					this->Edgemap.erase(itr3 ++);
					if((int)this->Resultmap.size() == this->k){
						break;
					}
				}
				if(this->Edgemap.size() == 0){
					return;
				}
				else{
					auto tmp_itr = this->Resultmap.end();
					tmp_itr --;
					while(itr3 != this->Edgemap.end()){
						if((*tmp_itr).first == (*itr3).first){
							this->Resultmap.insert(multimap<int,set<Edge*>>::value_type((*itr3).first, (*itr3).second));
							this->Edgemap.erase(itr3 ++);
						}
						else{
							break;
						}

					}
			    }
			}
		}
    }
#ifdef DEBUG
	cout << " out rm_old_result"  <<endl;
#endif
}

void basesolution::add_new_result(vector<int> &_L, vector<int>* _R){
	set<Edge*> tmp_biclique;
	int tmp_weight_cnt = 0;
	multimap< int,set<Edge*> >::iterator map_itr;
#ifdef DEBUG
	cout << "In add_new_result" <<endl;
	if((int)_R->size() < this->q){
		cerr << "_R.size() < this->q" <<endl;
		system("pause");
	}
	cout << "_R.size() = " <<(int)_R->size() <<endl;
	if((int)_L.size() != this->p){
		cerr << "_L.size() < this->p" <<endl;
		//getchar();
	}
	cout << "_L.size() = " <<(int)_L.size() <<endl;
#endif
	if((int)_R->size() >= this->q){
		vector< vector<int> > combs;//combinations of R
		if((int)_R->size() == this->q){
			combs.push_back(*(_R));
		}
		else{
			vector<int> comb(this->q);
			comb[this->q - 1] = this->snapshot[this->snapshot_size() - 1]->v;
			vector<int> rm_v_from_R;
			for(auto itr : *(_R)){
				if(itr != this->snapshot[this->snapshot_size() - 1]->v){
					rm_v_from_R.push_back(itr);
				}
			}
			this->find_combinations(combs, 0, 0, rm_v_from_R, comb);
		}
		//cerr << "combs.size() = " << combs.size() <<endl;
		for(int i = 0; i < (int)combs.size(); i ++){
			tmp_biclique.clear();
			tmp_weight_cnt = 0;
			//compute (p,q)-biclique and weight
			for(int j = 0; j < (int)_L.size(); j ++){
				cout << "_L[j] = " << _L[j] <<endl;
				cout << "this->one_hop_tb.size() = " << this->one_hop_tb.size() <<endl;
				for(auto itr : this->one_hop_tb[_L[j]]){
					sort(combs[i].begin(), combs[i].end());
						if(binary_search(combs[i].begin(), combs[i].end(), itr->v)){//binary_search of vector
							tmp_biclique.insert(itr);
							tmp_weight_cnt += itr->weight;
						}
						else{
							//cerr<< "compute _R or one_hop_tb error" <<endl;
							//system("pause");
						}
					}
			}

			if((int)this->Resultmap.size() < this->k){//如果top-k结果集长度小于k，此时edgemap必为空，则直接添加
				this->Resultmap.insert(multimap<int,set<Edge*>>::value_type(tmp_weight_cnt, tmp_biclique));
			}
			else{//top-k结果集长度大于等于k
				map_itr = this->Resultmap.end();
				map_itr--;
#ifdef DEBUG
				cout << "now this->Resultmap.size() = " << this->Resultmap.size() <<endl;
				cout << "tmp_weight_cnt = " << tmp_weight_cnt <<endl;
#endif
				if(tmp_weight_cnt < (*map_itr).first){//比较结果集中最后一个元素与候选biclique的大小
					this->Edgemap.insert(multimap<int,set<Edge*>>::value_type(tmp_weight_cnt, tmp_biclique));
				}
				else if(tmp_weight_cnt == (*map_itr).first){
					this->Resultmap.insert(multimap<int,set<Edge*>>::value_type(tmp_weight_cnt, tmp_biclique));
				}
				else{
					this->Resultmap.insert(multimap<int,set<Edge*>>::value_type(tmp_weight_cnt, tmp_biclique));
					int res_size = (int)this->Resultmap.size();
					int topk_loc = 0;
					for(map_itr = this->Resultmap.begin(); map_itr != this->Resultmap.end();){
						topk_loc ++;
						if(topk_loc != this->k){
							map_itr ++;
						}
						else{
							//auto tmp_itr = map_itr;
							if((*map_itr).first != (*++ map_itr).first){//如果topk中第k个元素与k+1个元素不相同
								while(map_itr != this->Resultmap.end()){
									this->Edgemap.insert(multimap<int,set<Edge*>>::value_type((*map_itr).first, (*map_itr).second));
									this->Resultmap.erase(map_itr ++);

								}
#ifdef DEBUG
								if((int)this->Resultmap.size() != this->k){
									cerr<< "return new result to edgemap failed " <<endl;
									system("pause");
								}
#endif
							}
						}

					}
				}
			}
		}
	}
	//cout << "OUT add_new_result" <<endl;
}

int basesolution::snapshot_size(){
	return this->snapshot.size();
}








