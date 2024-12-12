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

bool incrementVector(int size, vector<int>& currentVector) {
    int currentSize = currentVector.size(), j, maxVal;
    
    // Find the rightmost element that isn't at its maximum value
    for (j = currentSize - 1, maxVal = size - 1; j >= 0; j--) {
        if (currentVector[j] == maxVal) {
            maxVal--;
        } else {
            break;
        }
    }
    
    // If we're at the rightmost element, simply increment it
    if (j == currentSize - 1) {
        currentVector[j]++;
        return true;
    }
    
    // Otherwise, we need to reset subsequent elements
    int startMax = currentVector[j] + 1;
    while (j < currentSize) {
        currentVector[j++] = startMax++;
    }
    
    return true;
}

int computeTableLookup(int targetVertex, set<int>& vertexSet, vector<int>& allVertices) {
    int minimumCost = INT_MAX;
    
    // Special case for singleton set
    if (vertexSet.size() == 1) {
        auto singleVertex = *vertexSet.begin();
        set<pair<int,int>> temporaryPathSet;
        
        if (singleVertex == targetVertex) {
            C[{targetVertex, vertexSet}] = 0;
            Cset[{targetVertex, vertexSet}] = temporaryPathSet;
            return 0;
        } else {
            C[{targetVertex, vertexSet}] = d[targetVertex][singleVertex];
            temporaryPathSet.insert(dSet[{targetVertex, singleVertex}].begin(), 
                                    dSet[{targetVertex, singleVertex}].end());
            Cset[{targetVertex, vertexSet}] = temporaryPathSet;
            return d[targetVertex][singleVertex];
        }
    }
    
    // Convert set to vector for indexing
    vector<int> terminals(vertexSet.begin(), vertexSet.end());
    
    // First strategy: Choose a vertex from the set
    int firstChosenVertex = -1;
    int firstMinimumCost = INT_MAX;
    set<int> remainingSet1;
    
    for (auto intermediateVertex : vertexSet) {
        // Create a copy of the set without the intermediate vertex
        set<int> reducedSet(vertexSet);
        reducedSet.erase(intermediateVertex);
        
        // Compute cost recursively or retrieve from cache
        int intermediateCost = C.count({intermediateVertex, reducedSet}) ? 
                                C[{intermediateVertex, reducedSet}] : 
                                computeTableLookup(intermediateVertex, reducedSet, allVertices);
        
        // Update minimum cost path
        int potentialCost = intermediateCost + d[intermediateVertex][targetVertex];
        if (firstMinimumCost > potentialCost) {
            firstChosenVertex = intermediateVertex;
            firstMinimumCost = potentialCost;
        }
    }
    
    // Second strategy: Choose a vertex not in the set
    int secondChosenVertex = -1;
    int secondMinimumCost = INT_MAX;
    set<int> bestSubset;
    
    for (int candidateVertex : allVertices) {
        if (vertexSet.count(candidateVertex)) continue;
        
        // Generate all possible subset combinations
        vector<int> currentConfiguration(terminals.size(), -1);
        int totalCombinations = (1 << (terminals.size() - 1)) - 1;
        
        int localMinimumCost = INT_MAX;
        set<int> localBestSubset;
        
        for (int i = 0; i < totalCombinations; i++) {
            incrementVector(terminals.size(), currentConfiguration);
            
            // Construct subset based on configuration
            set<int> subset1, subset2;
            for (int j = currentConfiguration.size() - 1; 
                 j >= 0 && currentConfiguration[j] > -1; j--) {
                subset1.insert(terminals[currentConfiguration[j]]);
            }
            
            // Compute complement subset
            set_difference(vertexSet.begin(), vertexSet.end(), 
                           subset1.begin(), subset1.end(), 
                           inserter(subset2, subset2.begin()));
            
            // Compute total cost for this configuration
            int configurationCost = 
                (C.count({candidateVertex, subset1}) ? C[{candidateVertex, subset1}] : 
                 computeTableLookup(candidateVertex, subset1, allVertices)) +
                (C.count({candidateVertex, subset2}) ? C[{candidateVertex, subset2}] : 
                 computeTableLookup(candidateVertex, subset2, allVertices));
            
            // Update local minimum
            if (configurationCost < localMinimumCost) {
                localMinimumCost = configurationCost;
                localBestSubset = subset1;
            }
        }
        
        // Update global minimum
        int potentialCost = localMinimumCost + d[candidateVertex][targetVertex];
        if (secondMinimumCost > potentialCost) {
            secondMinimumCost = potentialCost;
            secondChosenVertex = candidateVertex;
            bestSubset = localBestSubset;
        }
    }
    
    // Choose the better strategy
    set<pair<int, int>> finalPathSet;
    set<int> remainingSet;
    
    if (firstMinimumCost < secondMinimumCost) {
        minimumCost = firstMinimumCost;
        remainingSet = vertexSet;
        remainingSet.erase(firstChosenVertex);
        
        // Merge path sets
        finalPathSet.insert(Cset[{firstChosenVertex, remainingSet}].begin(), 
                            Cset[{firstChosenVertex, remainingSet}].end());
        finalPathSet.insert(dSet[{firstChosenVertex, targetVertex}].begin(), 
                            dSet[{firstChosenVertex, targetVertex}].end());
    } else {
        minimumCost = secondMinimumCost;
        
        // Compute the complement subset
        set<int> complementSet;
        set_difference(vertexSet.begin(), vertexSet.end(), 
                       bestSubset.begin(), bestSubset.end(), 
                       inserter(complementSet, complementSet.begin()));
        
        // Merge path sets from both subsets
        finalPathSet.insert(Cset[{secondChosenVertex, bestSubset}].begin(), 
                            Cset[{secondChosenVertex, bestSubset}].end());
        finalPathSet.insert(Cset[{secondChosenVertex, complementSet}].begin(), 
                            Cset[{secondChosenVertex, complementSet}].end());
        finalPathSet.insert(dSet[{secondChosenVertex, targetVertex}].begin(), 
                            dSet[{secondChosenVertex, targetVertex}].end());
    }
    
    // Update global caches
    C[{targetVertex, vertexSet}] = minimumCost;
    Cset[{targetVertex, vertexSet}] = finalPathSet;
    
    return minimumCost;
}

int minVal = INT_MAX;
set<pair<int,int>> solPSet;
	
void parseGraphSection(vector<vector<Edge>>& graph, map<pair<int,int>, int>& W, 
                       vector<long>& incident, set<long>& weightSet) {
    string dummy;
    long m, n, u, v, w;
    
    cin >> dummy >> n;
    cin >> dummy >> m;
    
    incident.resize(n+1, 0);
    graph.resize(n+1);
    
    for (long i = 0; i < m; i++) {
        cin >> dummy >> u >> v >> w;
        
        weightSet.insert(w);
        incident[u] += 1;
        incident[v] += 1;
        
        graph[u].push_back(Edge(v, w));
        graph[v].push_back(Edge(u, w));
        W[make_pair(u, v)] = w;
        W[make_pair(v, u)] = w;
    }
    cin >> dummy;
}

void parseTerminalsSection(vector<int>& terminals) {
    string dummy;
    long t, u;
    
    cin >> dummy >> t;
    for (long i = 0; i < t; i++) {
        cin >> dummy >> u;
        terminals.push_back(u);
    }
    cin >> dummy;
}

void parseTreeSection() {
    string dummy, line;
    long b, val;
    
    cin >> dummy >> dummy >> dummy;
    cin >> b;
    
    cin >> dummy >> dummy >> ws;
    
    for (long i = 0; i < b; i++) {
        getline(cin, line);
        stringstream sstream(line);
        if (sstream >> dummy, dummy == "b") {
            while (sstream >> val) {
                // Process values if needed
            }
        }
    }
    
    long tu, tv;
    for (long i = 0; i < b - 1; i++) {
        cin >> tu >> tv;
    }
    cin >> dummy; // end
}

void findMinimumSteinerTree(vector<int>& terminals, vector<vector<Edge>>& graph, 
                             map<pair<int,int>, int>& W) {
    int n = graph.size();
    
    // Compute all-pairs shortest paths
    d = dijkstra(graph, W);
    
    set<int> terminalSet(terminals.begin(), terminals.end());
    vector<int> V;
    for (int i = 1; i < n; i++)
        V.push_back(i);
    
    minVal = INT_MAX;
    solPSet.clear();
    
    for (auto v : terminals) {
        terminalSet.erase(v);
        int val = computeTableLookup(v, terminalSet, V);
        
        if (val < minVal) {
            minVal = val;
            solPSet.clear();
            solPSet.insert(Cset[{v, terminalSet}].begin(), Cset[{v, terminalSet}].end());
        }
        terminalSet.insert(v);
    }
}

void outputResults(map<pair<int,int>, int>& W) {
    cout << "VALUE " << minVal << endl;
    int totalWeight = 0;
    for (auto u : solPSet) {
        cout << u.first << " " << u.second << endl;
        totalWeight += W[{u.first, u.second}];
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    
    vector<vector<Edge>> graph;
    map<pair<int,int>, int> W;
    vector<int> terminals;
    
    set<long> weightSet;
    vector<long> incident;
    
    string code, type;
    
    while (cin >> code >> type) {
        if (code == "SECTION") {
            if (type == "Graph") {
                parseGraphSection(graph, W, incident, weightSet);
            }
            else if (type == "Terminals") {
                parseTerminalsSection(terminals);
            }
            else if (type == "Tree") {
                parseTreeSection();
            }
        }
    }
    
    findMinimumSteinerTree(terminals, graph, W);
    outputResults(W);
    
    return 0;
}