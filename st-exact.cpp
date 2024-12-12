#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <sstream>
#include <algorithm>
#include <climits>
#include <signal.h>
#include <unistd.h>
#include <cstring>
#include <cmath>
#include <omp.h>
#include <queue>
#include <iterator>

using namespace std;
class Edge { 
	
public:
	int to;
	int length; 
	
	Edge(){}
	~Edge(){}
	Edge(int t, int l){
		to = t; length = l;
	}
	bool operator < (const Edge& e){
		return length < e.length;
	}
};


// recordes the edges on path(u,v)
map<pair<int,int> , set<pair<int,int>> > dSet; 

vector<vector<int>> dijkstra(
    vector<vector<Edge>> &graph,
    map<pair<int, int>, int> W) {
  int n = graph.size(); // Number of vertices in the graph
  vector<int> empty(n, INT_MAX / 2 - 3); // Ensures sum of two INT_MAX stays within limits
  vector<vector<int>> dist(n, empty);

  // Initialize distance for direct edges based on weights in W
  for (auto pairW : W) {
    int u = pairW.first.first;
    int v = pairW.first.second;
    int w = pairW.second;
    dist[u][v] = w;
    dist[v][u] = w;

    pair<int, int> tmpPair;
    if (u < v) // To ensure unique edge pairs
      tmpPair = {u, v};
    else
      tmpPair = {v, u};

    dSet[{u, v}].insert(tmpPair); // Record edges in dSet
    dSet[{v, u}].insert(tmpPair);
  }

  // Perform Dijkstra's algorithm for each vertex as the source
  for (int src = 1; src < n; src++) {
    vector<int> minDist(n, INT_MAX / 2 - 3);
    minDist[src] = 0;
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    pq.push({0, src});

    map<int, set<pair<int, int>>> localDSet;

    while (!pq.empty()) {
      int currDist = pq.top().first;
      int u = pq.top().second;
      pq.pop();

      if (currDist > minDist[u])
        continue;

      for (const Edge &edge : graph[u]) {
        int v = edge.to;
        int weight = edge.length;

        if (minDist[u] + weight < minDist[v]) {
          minDist[v] = minDist[u] + weight;
          pq.push({minDist[v], v});

          // Update the path information
          localDSet[v].clear();
          localDSet[v].insert(localDSet[u].begin(), localDSet[u].end());
          localDSet[v].insert({u, v});
        }
      }
    }

    // Update the global distance and dSet
    for (int dest = 1; dest < n; dest++) {
      dist[src][dest] = minDist[dest];
      dist[dest][src] = minDist[dest];

      dSet[{src, dest}].clear();
      dSet[{src, dest}].insert(localDSet[dest].begin(), localDSet[dest].end());

      dSet[{dest, src}].clear();
      dSet[{dest, src}].insert(localDSet[dest].begin(), localDSet[dest].end());
    }
  }

  return dist;
}


map<pair<int, set<int>>,int> C; 
map<pair<int, set<int>>,set<pair<int,int>> > Cset; 

vector <vector<int>> d;

bool incVec(int size, vector<int>& curV){
	int cSize = curV.size(), j, max;
	for(j = cSize-1, max = size-1; j >= 0; j--) {
		if(curV[j] == max) {
			max--;
		}
		else {
			break;
		}
	}
	if(j == cSize-1) {
		curV[j]++;
	}
	else {
		while (j > 0 && curV[j] >= max) {
			j--;
			max--;
		}
		max = curV[j]+1;
		while(j < cSize) {
			curV[j++] = max++;
		}
	}
	return true;
}

int computeTab(int v, set<int> & X, vector <int> & V){
	
	int minVal =INT_MAX;
	if(X.size() ==1){
		auto x0 = *(X.begin());
		set<pair<int,int>> tmpPSet; 
		if(x0 == v){
			C.insert({{v,X}, 0});
			Cset.insert({{v,X}, tmpPSet});
			return 0;
		}	
		else{ // x0 is some u
			C.insert({{v,X}, d[v][x0]});
			tmpPSet.insert(dSet[{v,x0}].begin(), dSet[{v,x0}].end());
			
			Cset.insert({{v,X}, tmpPSet});
		}
		return d[v][x0];
	}
	
	vector<int> terminals(X.begin(), X.end()); // making set x as vector!
	
	set<int> X1, X2;
	X1.insert(X.begin(), X.end());
	int v1 = -1;
	int minFirst = INT_MAX;
	for(auto u: X){
		X1.erase(u);
		int cost = -1;
		if(C.find({u,X1}) == C.end()) {
			cost = computeTab(u, X1, V);
		}
		else {
			cost = C[{u,X1}];
		}
		if(minFirst > cost+d[u][v]) {
			v1 = u;
			minFirst = cost + d[u][v];
		}
		X1.insert(u);
	}
	
	int v2 = -1;
	int minSecond = INT_MAX;
	set<int> Xprime;
	int it = -1;
	for(auto u: V) {
		
		if(X.find(u) != X.end()) continue;
		vector<int> curV(terminals.size(),-1); 
		int size = (1 << (terminals.size() - 1)) -1 ; // 2^(n-1) -1
		
		set<int> Xlocal;
		int minLocal = INT_MAX;
		for(int i= 0; i < size; i++) {
			incVec(terminals.size(), curV);
			X1.clear();
			for(int j = curV.size()-1; j >= 0 && curV[j] > -1; j--) {
				X1.insert(terminals[curV[j]]);
			}
			X2.clear();
			set_difference(X.begin(), X.end(), X1.begin(), X1.end(), inserter(X2, X2.begin()));
			int cost = 0;
			if(C.find({u,X1}) == C.end()) {
				cost = computeTab(u, X1, V);
			}
			else {
				cost = C[{u,X1}];
			}
			if(C.find({u,X2}) == C.end()) {
				cost += computeTab(u, X2, V);
			}
			else {
				cost += C[{u,X2}];
			}
			if(cost < minLocal) {
				minLocal = cost;
				Xlocal.clear();
				Xlocal.insert(X1.begin(), X1.end());
			}
 		}
 		if(minSecond > minLocal + d[u][v]) {
			minSecond = minLocal + d[u][v];
			v2 = u;
			Xprime.clear();
			Xprime.insert(Xlocal.begin(), Xlocal.end());
		}
	}
	set<pair<int, int>> P;
	set<int> Xlocal(X.begin(), X.end());
	if(minFirst < minSecond) {
		minVal = minFirst;
		Xlocal.erase(v1);
		auto Pset = Cset[{v1,Xlocal}];
		for(auto e: Pset) {
			P.insert(e);
		}
		for(auto e: dSet[{v1,v}]) {
			P.insert(e);
		}
	}
	else {
		minVal = minSecond;
		Xlocal.clear();
		set_difference(X.begin(), X.end(), Xprime.begin(), Xprime.end(), inserter(Xlocal, Xlocal.begin()));		
		auto Pset = Cset[{v2,Xprime}];
		for(auto e: Pset) {
			P.insert(e);
		}
		Pset = Cset[{v2,Xlocal}];
		for(auto e: Pset) {
			P.insert(e);
		}
		for(auto e: dSet[{v2,v}]) {
			P.insert(e);
		}
	}
	C[{v,X}] = minVal;
	Cset.insert({{v,X}, P});	
	return minVal;
}

int minVal = INT_MAX;
set<pair<int,int>> solPSet;

int main(){
	
 
    ios_base::sync_with_stdio(false);

	vector< vector<Edge> > graph;
	map<pair<int,int> , int> W;
	vector <int> terminals;
	string code, type, dummy;
	

	set<long> weightSet;
	vector<long> incident;
	 
	while( cin>> code >> type ){

		if(code == "SECTION" && type =="Graph"){
			long m, n;
			long u, v, w;
			cin >> dummy >> n;
			cin >> dummy >> m;
			
			incident.resize(n+1, 0) ;
			
			graph.resize(n+1); // coz graph has from index 0. where as challege its 1
			for(long i=0; i < m; i++){
				cin>> dummy >> u >> v >> w;
				
				weightSet.insert(w);
				incident[u]+=1;
				incident[v]+=1;
				
				graph[u].push_back(Edge(v,w));
				graph[v].push_back(Edge(u,w));
				W[make_pair(u,v)]=w;
				W[make_pair(v,u)]=w;
			}
			cin >> dummy;
		}
		else if(code == "SECTION" && type =="Terminals"){
			long t, u;
			cin >> dummy >> t;
			for(long i=0; i < t; i++){
				cin>> dummy >> u;
				//cout << "T" << u << endl;
				terminals.push_back(u);
			}
			cin >> dummy;
		}
		else if(code == "SECTION" && type =="Tree"){
			
			cin >> dummy >> dummy >> dummy;
			long b, val ; cin >> b; 
			
			cin >> dummy >> dummy >> ws;
			
			for(long i=0; i < b; i++){
				string line;
				getline(cin, line); stringstream sstream(line);
				if(sstream >> dummy, dummy=="b"){
					while(sstream >> val){
						//cout << val << " " ;
					}
				}
			}
			long tu, tv;
			for(long i=0; i < b-1; i++){ // b-1 edges is Td
				cin >>  tu >> tv;
			}
			cin >> dummy; // end
		}
	}
	
	d = dijkstra(graph, W); //n^2logn, n^2 space
	int n = graph.size();

	set <int> terminalSet(terminals.begin(), terminals.end());
	vector <int> V;
	for(int i=1; i < n; i++)
		V.push_back(i);

	for(auto v : terminals){
		terminalSet.erase(v);
		int val = computeTab(v, terminalSet, V);

		if(val < minVal){
			minVal = val;
			solPSet.clear();
			solPSet.insert(Cset[{v,terminalSet}].begin(), Cset[{v,terminalSet}].end());
		}
		terminalSet.insert(v);
	}
	cout << "VALUE " << minVal << endl;
	int y = 0;
	for(auto u: solPSet) {
		cout<<u.first<<" "<<u.second<< endl;
		y += W[{u.first,u.second}];
	}
	
	return 0;
}
	
