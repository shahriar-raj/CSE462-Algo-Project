#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <climits>
#include <algorithm>
#include <string>
#include <sstream>
#include <random>

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
vector<int> dijkstraShortestPath(const vector<vector<Edge>>& graph, int start) {
    int n = (int)graph.size();
    vector<int> dist(n, INT_MAX);
    vector<int> prev(n, -1);

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

long long buildSteinerTreeFromStart(int startVertex, const vector<int>& terminals,
                                    const vector<vector<Edge>>& graph,
                                    const map<pair<int,int>, int>& originalW,
                                    map<pair<int,int>, int>& bestEdges) {
    // Copy W so we don't modify the original map
    map<pair<int,int>, int> W = originalW;

    set<int> steinerTreeVertices;
    steinerTreeVertices.insert(startVertex);

    map<pair<int,int>, int> steinerTreeEdges;
    long long totalWeight = 0;

    vector<int> remainingTerminals = terminals;
    remainingTerminals.erase(
        remove(remainingTerminals.begin(), remainingTerminals.end(), startVertex),
        remainingTerminals.end()
    );

    while (!remainingTerminals.empty()) {
        int nearestTerminal = -1;
        int minPathLength = INT_MAX;
        int connectingVertex = -1;

        // Temporary best edges for this iteration
        map<pair<int,int>, int> iterationBestEdges;

        for (int treeVertex : steinerTreeVertices) {
            vector<int> prev = dijkstraShortestPath(graph, treeVertex);

            for (int terminal : remainingTerminals) {
                // Reconstruct path
                int current = terminal;
                int pathLength = 0;
                vector<pair<int,int>> pathEdges;

                while (current != -1 && current != treeVertex) {
                    int parent = prev[current];
                    if (parent == -1) break;

                    int a = min(parent, current);
                    int b = max(parent, current);
                    auto it = W.find({a,b});
                    if (it == W.end()) {
                        // edge not found, break
                        pathEdges.clear();
                        break;
                    }
                    int edgeWeight = it->second;

                    pathLength += edgeWeight;
                    pathEdges.push_back({a,b});
                    current = parent;
                }

                if (current == treeVertex && !pathEdges.empty() && pathLength < minPathLength) {
                    minPathLength = pathLength;
                    nearestTerminal = terminal;
                    connectingVertex = treeVertex;

                    iterationBestEdges.clear();
                    for (auto &ed : pathEdges) {
                        iterationBestEdges[ed] = W.at(ed);
                    }
                }
            }
        }

        if (nearestTerminal == -1) break;

        for (auto& edge : iterationBestEdges) {
            steinerTreeVertices.insert(edge.first.first);
            steinerTreeVertices.insert(edge.first.second);
            totalWeight += edge.second;
            steinerTreeEdges[edge.first] = edge.second;
        }

        remainingTerminals.erase(
            remove(remainingTerminals.begin(), remainingTerminals.end(), nearestTerminal),
            remainingTerminals.end()
        );
    }

    // Return result
    if (totalWeight < LLONG_MAX) {
        bestEdges = steinerTreeEdges;
    }
    return totalWeight;
}


void findMinimumSteinerTree(const vector<int>& terminals,
                            const vector<vector<Edge>>& graph,
                            map<pair<int,int>, int>& W) {
    if (terminals.empty()) return;

    // We will try multiple attempts with different random starting terminals
    // and pick the best solution.
    static random_device rd;
    static mt19937 gen(rd());
    // If you have at least one terminal, pick among them randomly.
    if ((int)terminals.size() == 1) {
        // Only one terminal, must start there.
        map<pair<int,int>, int> bestEdges;
        long long cost = buildSteinerTreeFromStart(terminals[0], terminals, graph, W, bestEdges);
        minVal = cost;
        W = bestEdges;
        return;
    }

    int num_attempts = 5; // number of random attempts
    long long best_cost = LLONG_MAX;
    map<pair<int,int>, int> best_solution;

    uniform_int_distribution<int> dist(0,(int)terminals.size()-1);

    for (int attempt = 0; attempt < num_attempts; attempt++) {
        int startTerminal = terminals[dist(gen)];
        map<pair<int,int>, int> attempt_edges;
        long long cost = buildSteinerTreeFromStart(startTerminal, terminals, graph, W, attempt_edges);
        if (cost < best_cost) {
            best_cost = cost;
            best_solution = attempt_edges;
        }
    }

    minVal = best_cost;
    W = best_solution;
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
        W[make_pair(min(u,v), max(u,v))] = w;
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
