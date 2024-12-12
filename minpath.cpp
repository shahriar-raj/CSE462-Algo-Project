#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <climits>
#include <algorithm>
#include <string>
#include <sstream>

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

// Global variable to store the minimum tree value
long long minVal = LLONG_MAX;

// Dijkstra's shortest path algorithm
vector<int> dijkstraShortestPath(const vector<vector<Edge>>& graph, int start, const vector<int>& terminals) {
    int n = graph.size();
    vector<int> dist(n, INT_MAX);
    vector<int> prev(n, -1);
    vector<bool> inTerminals(n, false);
    
    for (int terminal : terminals) {
        inTerminals[terminal] = true;
    }
    
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    dist[start] = 0;
    pq.push({0, start});
    
    while (!pq.empty()) {
        int u = pq.top().second;
        int d = pq.top().first;
        pq.pop();
        
        if (d > dist[u]) continue;
        
        for (const Edge& edge : graph[u]) {
            int v = edge.to;
            int newDist = d + edge.length;
            
            if (newDist < dist[v]) {
                dist[v] = newDist;
                prev[v] = u;
                pq.push({newDist, v});
            }
        }
    }
    
    return prev;
}

void findMinimumSteinerTree(const vector<int>& terminals, 
                             const vector<vector<Edge>>& graph, 
                             map<pair<int,int>, int>& W) {
    if (terminals.empty()) return;
    
    // Start with an arbitrary terminal
    int startVertex = terminals[0];
    
    // Set to keep track of vertices in the Steiner tree
    set<int> steinerTreeVertices;
    steinerTreeVertices.insert(startVertex);
    
    // Map to store the edges of the Steiner tree
    map<pair<int,int>, int> steinerTreeEdges;
    long long totalWeight = 0;
    
    // Copy of terminals to track which are not yet in the tree
    vector<int> remainingTerminals = terminals;
    remainingTerminals.erase(
        remove(remainingTerminals.begin(), remainingTerminals.end(), startVertex), 
        remainingTerminals.end()
    );
    
    // Iteratively add the nearest terminal to the tree
    while (!remainingTerminals.empty()) {
        int nearestTerminal = -1;
        int minPathLength = INT_MAX;
        int connectingVertex = -1;
        
        // Find the shortest path from current tree to a remaining terminal
        for (int treeVertex : steinerTreeVertices) {
            vector<int> prev = dijkstraShortestPath(graph, treeVertex, remainingTerminals);
            
            for (int terminal : remainingTerminals) {
                // Reconstruct path
                int current = terminal;
                int pathLength = 0;
                vector<pair<int,int>> pathEdges;
                
                while (current != -1 && current != treeVertex) {
                    int parent = prev[current];
                    if (parent == -1) break;
                    
                    // Find the weight of this edge
                    int edgeWeight = 0;
                    auto it = W.find({parent, current});
                    if (it != W.end()) {
                        edgeWeight = it->second;
                    }
                    
                    pathLength += edgeWeight;
                    pathEdges.push_back({parent, current});
                    current = parent;
                }
                
                // Update if this is the shortest path to a terminal
                if (current == treeVertex && pathLength < minPathLength) {
                    minPathLength = pathLength;
                    nearestTerminal = terminal;
                    connectingVertex = treeVertex;
                    
                    // Store the path edges
                    steinerTreeEdges.clear();
                    for (auto edge : pathEdges) {
                        steinerTreeEdges[edge] = W.at(edge);
                    }
                }
            }
        }
        
        // If no path found, break
        if (nearestTerminal == -1) break;
        
        // Add the path to the Steiner tree
        for (auto& edge : steinerTreeEdges) {
            steinerTreeVertices.insert(edge.first.first);
            steinerTreeVertices.insert(edge.first.second);
            totalWeight += edge.second;
        }
        
        // Remove the added terminal
        remainingTerminals.erase(
            remove(remainingTerminals.begin(), remainingTerminals.end(), nearestTerminal), 
            remainingTerminals.end()
        );
    }
    
    // Update global minimum value
    minVal = totalWeight;
    
    // Update W with the Steiner tree edges
    W = steinerTreeEdges;
}

// Existing parsing functions (kept as in the original code)
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
 
void outputResults(map<pair<int,int>, int>& W) { 
    cout << "VALUE " << minVal << endl;
    
    // Print out edges of the Steiner tree
    cout << "EDGES" << endl;
    for (const auto& edge : W) {
        cout << edge.first.first << " " << edge.first.second << endl;
    }
    cout << "END" << endl;
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