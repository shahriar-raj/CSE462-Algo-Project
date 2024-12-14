#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <climits>
#include <algorithm>
#include <string>
#include <sstream>
#include <unordered_set>

using namespace std;


// Struct Definitions


struct Graph {
    // Adjacency list: node -> [(neighbor, weight)]
    vector<vector<pair<int, int>>> adj;
    int num_nodes;

    Graph(int n) : adj(n + 1), num_nodes(n) {}

    // Add undirected edge
    void add_edge(int u, int v, int weight) {
        adj[u].emplace_back(v, weight);
        adj[v].emplace_back(u, weight);
    }
};


// Global Variables


static long long minVal = LLONG_MAX;
static vector<int> global_terminals;
static map<pair<int, int>, int> bestSolutionEdges;
static vector<bool> is_terminal;
static int n; // Number of nodes

// The global edge map W stores all edges of the input graph in normalized form:
// For an edge (u,v), we always store it as (min(u,v), max(u,v)).
static map<pair<int, int>, int> original_W;


// Function Prototypes


// Dijkstra's shortest path: returns prev[] for path reconstruction
vector<int> dijkstraShortestPath(const Graph& graph, int start);

// Check if all terminals are connected in a given set of edges
bool terminalsConnected(const set<pair<int, int>>& edges, const vector<int>& terminals);

// Compute cost of a solution given its edges and original_W
long long computeCost(const set<pair<int, int>>& edges);

// Build a solution starting from a given terminal
map<pair<int, int>, int> buildSolutionFromTerminal(int startVertex,
                                                   const vector<int>& terminals,
                                                   const Graph& graph);

// Enhanced Local Search: Remove and Add non-terminal nodes to reduce cost
bool local_search_improved(map<pair<int, int>, int>& W_map, const vector<int>& terminals, const Graph& graph);

// Parse different sections of the input
void parseGraphSection(Graph& graph, map<pair<int, int>, int>& W,
                       vector<long>& incident, set<long>& weightSet);

void parseTerminalsSection(vector<int>& terminals);

void parseTreeSection();

// Output the results
void outputResults(const map<pair<int, int>, int>& W_map);


// Function Implementations


// Dijkstra's shortest path: returns prev[] for path reconstruction
vector<int> dijkstraShortestPath(const Graph& graph, int start) {
    int N = graph.num_nodes;
    vector<int> dist(N + 1, INT_MAX);
    vector<int> prev(N + 1, -1);

    dist[start] = 0;
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    pq.push({0, start});

    while (!pq.empty()) {
        int u = pq.top().second;
        int d = pq.top().first;
        pq.pop();

        if (d > dist[u]) continue;

        for (const auto& edge : graph.adj[u]) {
            int v = edge.first;
            int weight = edge.second;
            int newDist = d + weight;

            if (newDist < dist[v]) {
                dist[v] = newDist;
                prev[v] = u;
                pq.push({newDist, v});
            }
        }
    }

    return prev;
}

// Check if all terminals are connected in a given set of edges
bool terminalsConnected(const set<pair<int, int>>& edges, const vector<int>& terminals) {
    if (terminals.empty()) return true;

    // Build adjacency list from edges
    vector<vector<int>> adj(n + 1);
    for (const auto& ed : edges) {
        int u = ed.first, v = ed.second;
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    // BFS to check connectivity
    int start = terminals[0];
    unordered_set<int> visited;
    queue<int> q;
    q.push(start);
    visited.insert(start);

    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (const auto& neighbor : adj[u]) {
            if (!visited.count(neighbor)) {
                visited.insert(neighbor);
                q.push(neighbor);
            }
        }
    }

    // Check if all terminals are visited
    for (const auto& term : terminals) {
        if (!visited.count(term)) return false;
    }
    return true;
}

// Compute cost of a solution given its edges and original_W
long long computeCost(const set<pair<int, int>>& edges) {
    long long cost = 0;
    for (const auto& e : edges) {
        // e is already normalized
        cost += original_W.at(e);
    }
    return cost;
}

// Build a solution starting from a given terminal
map<pair<int, int>, int> buildSolutionFromTerminal(int startVertex,
                                                   const vector<int>& terminals,
                                                   const Graph& graph) {
    set<int> steinerTreeVertices;
    steinerTreeVertices.insert(startVertex);

    map<pair<int, int>, int> steinerTreeEdges;
    vector<int> remainingTerminals = terminals;
    // Remove startVertex from remainingTerminals if present
    remainingTerminals.erase(
        remove(remainingTerminals.begin(), remainingTerminals.end(), startVertex),
        remainingTerminals.end()
    );

    while (!remainingTerminals.empty()) {
        int nearestTerminal = -1;
        int minPathLength = INT_MAX;

        map<pair<int, int>, int> bestPathEdges;

        for (const int treeVertex : steinerTreeVertices) {
            vector<int> prev = dijkstraShortestPath(graph, treeVertex);

            for (const int terminal : remainingTerminals) {
                int current = terminal;
                int pathLength = 0;
                vector<pair<int, int>> pathEdges;

                // Reconstruct path from terminal to treeVertex
                while (current != -1 && current != treeVertex) {
                    int parent = prev[current];
                    if (parent == -1) break;

                    int a = min(parent, current);
                    int b = max(parent, current);
                    pair<int, int> ed(a, b);

                    // Ensure edge is in original_W
                    if (original_W.find(ed) == original_W.end()) {
                        // Edge not found, invalid path
                        pathEdges.clear();
                        break;
                    }

                    pathLength += original_W.at(ed);
                    pathEdges.emplace_back(ed);
                    current = parent;
                }

                if (current == treeVertex && !pathEdges.empty() && pathLength < minPathLength) {
                    minPathLength = pathLength;
                    nearestTerminal = terminal;
                    bestPathEdges.clear();
                    for (const auto& edge : pathEdges) {
                        bestPathEdges[edge] = original_W.at(edge);
                    }
                }
            }
        }

        if (nearestTerminal == -1) {
            // No path found, cannot connect all terminals
            break;
        }

        // Add the best path edges to the Steiner Tree
        for (const auto& edge : bestPathEdges) {
            steinerTreeEdges[edge.first] = edge.second;
            steinerTreeVertices.insert(edge.first.first);
            steinerTreeVertices.insert(edge.first.second);
        }

        // Remove the connected terminal from remainingTerminals
        remainingTerminals.erase(
            remove(remainingTerminals.begin(), remainingTerminals.end(), nearestTerminal),
            remainingTerminals.end()
        );
    }

    return steinerTreeEdges;
}

// Local Search: Remove and Add non-terminal nodes to reduce cost
bool local_search_improved(map<pair<int, int>, int>& W_map, const vector<int>& terminals, const Graph& graph) {
    set<pair<int, int>> edgeSet;
    for (const auto& e : W_map) edgeSet.emplace(e.first);

    long long best_cost = computeCost(edgeSet);

    // Identify nodes in solution
    unordered_set<int> nodes_in_solution;
    for (const auto& ed : edgeSet) {
        nodes_in_solution.insert(ed.first);
        nodes_in_solution.insert(ed.second);
    }

    // **Phase 1: Attempt to Remove Non-Terminal Nodes**
    for (const int node : nodes_in_solution) {
        if (is_terminal[node]) continue; // Skip terminals

        // Attempt to remove this node by removing all connected edges
        set<pair<int, int>> new_solution = edgeSet;
        vector<pair<int, int>> to_remove;
        for (const auto& ed : edgeSet) {
            if (ed.first == node || ed.second == node) {
                to_remove.emplace_back(ed);
            }
        }
        for (const auto& r : to_remove) {
            new_solution.erase(r);
        }

        // Check if terminals remain connected
        if (terminalsConnected(new_solution, terminals)) {
            long long new_cost = computeCost(new_solution);
            if (new_cost < best_cost) {
                // Found improvement by removal
                map<pair<int, int>, int> newW;
                for (const auto& ed : new_solution) {
                    newW[ed] = original_W.at(ed);
                }
                W_map = newW;
                minVal = new_cost;
                bestSolutionEdges = newW;
                return true; // Early exit on improvement
            }
        }
    }

    // **Phase 2: Attempt to Add Non-Terminal Nodes**
    for (int node = 1; node <= graph.num_nodes; node++) {
        if (nodes_in_solution.find(node) != nodes_in_solution.end() || is_terminal[node]) {
            continue; // Skip already included nodes and terminals
        }

        // Identify connections from this node to the current solution
        vector<pair<int, int>> candidate_edges;
        for (const auto& edge : graph.adj[node]) {
            if (nodes_in_solution.find(edge.first) != nodes_in_solution.end()) {
                int a = min(node, edge.first);
                int b = max(node, edge.first);
                candidate_edges.emplace_back(a, b);
            }
        }

        // Heuristic: Add node if it connects to at least two nodes in the solution
        if (candidate_edges.size() >= 2) { // Threshold can be adjusted
            set<pair<int, int>> new_solution = edgeSet;
            for (const auto& ed : candidate_edges) {
                new_solution.emplace(ed);
            }

            // Check if terminals remain connected
            if (terminalsConnected(new_solution, terminals)) {
                long long new_cost = computeCost(new_solution);
                if (new_cost < best_cost) {
                    // Found improvement by addition
                    map<pair<int, int>, int> newW;
                    for (const auto& ed : new_solution) {
                        newW[ed] = original_W.at(ed);
                    }
                    W_map = newW;
                    minVal = new_cost;
                    bestSolutionEdges = newW;
                    return true; // Early exit on improvement
                }
            }
        }
    }

    return false; // No improvement found
}

// Parse the Graph section from input
void parseGraphSection(Graph& graph, map<pair<int, int>, int>& W,
                       vector<long>& incident, set<long>& weightSet) {
    string dummy;
    long m, u, v, w;
    cin >> dummy >> n; // Read "Nodes <n>"
    cin >> dummy >> m; // Read "Edges <m>"

    // Resize the adjacency list to accommodate all nodes
    graph.adj.assign(n + 1, vector<pair<int, int>>());
    graph.num_nodes = n;

    incident.resize(n + 1, 0); // Ensure 'incident' is also correctly sized

    for (long i = 0; i < m; i++) {
        cin >> dummy >> u >> v >> w; // Read "Edge u v w"

        // Validate node indices
        if (u < 1 || u > n || v < 1 || v > n) {
            cerr << "Error: Node indices out of bounds. Nodes should be between 1 and " << n << ".\n";
            exit(1);
        }

        weightSet.emplace(w);
        incident[u] += 1;
        incident[v] += 1;

        graph.add_edge(u, v, w);

        pair<int, int> ed(min(u, v), max(u, v));
        W[ed] = w;
    }
    cin >> dummy; // Read "END"
}

// Parse the Terminals section from input
void parseTerminalsSection(vector<int>& terminals) {
    string dummy;
    long t, u;

    cin >> dummy >> t; // Read "Terminals <t>"
    for (long i = 0; i < t; i++) {
        cin >> dummy >> u; // Read "Terminal <u>"
        terminals.emplace_back(u);
    }
    cin >> dummy; // Read "END"
}

// Parse the Tree section from input (if needed)
void parseTreeSection() {
    string dummy, line;
    long b, val;

    cin >> dummy >> dummy >> dummy; // Read "Tree SECTION <...>"
    cin >> b; // Read number of edges or some value

    cin >> dummy >> dummy >> ws; // Read more data as needed

    for (long i = 0; i < b; i++) {
        getline(cin, line);
        // Process if needed
    }

    long tu, tv;
    for (long i = 0; i < b - 1; i++) {
        cin >> tu >> tv;
    }
    cin >> dummy; // Read "END"
}

// Output the results
void outputResults(const map<pair<int, int>, int>& W_map) {
    cout << "VALUE " << minVal << endl;

    cout << "EDGES" << endl;
    for (const auto& edge : W_map) {
        cout << edge.first.first << " " << edge.first.second << endl;
    }
    cout << "END" << endl;
}


// Main Function


int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL); // For faster input

    // Variables to store graph information
    Graph graph(0); // Temporary initialization; will set actual number after parsing
    map<pair<int, int>, int> W;
    vector<int> terminals;

    set<long> weightSet;
    vector<long> incident;

    string code, type;

    // Parse the input
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

    // Check if terminals are provided
    if (terminals.empty()) {
        cerr << "Error: No terminals provided.\n";
        exit(1);
    }

    // Initialize global variables
    original_W = W; // Keep a global reference to all edges
    global_terminals = terminals;
    is_terminal.assign(n + 1, false);
    for (const auto& term : terminals) {
        if (term < 1 || term > n) {
            cerr << "Error: Terminal node " << term << " is out of bounds.\n";
            exit(1);
        }
        is_terminal[term] = true;
    }

    // Try multiple initial solutions from different terminals
    long long best_init_cost = LLONG_MAX;
    map<pair<int, int>, int> best_init_solution;

    for (const auto& startT : terminals) {
        auto candidate = buildSolutionFromTerminal(startT, terminals, graph);
        set<pair<int, int>> candEdges;
        for (const auto& ed : candidate) candEdges.emplace(ed.first);
        if (terminalsConnected(candEdges, terminals)) {
            long long cost = computeCost(candEdges);
            if (cost < best_init_cost) {
                best_init_cost = cost;
                best_init_solution = candidate;
            }
        }
    }

    if (best_init_solution.empty()) {
        // Fallback if no solution found
        best_init_solution = buildSolutionFromTerminal(terminals[0], terminals, graph);
        set<pair<int, int>> edgeSet;
        for (const auto& ed : best_init_solution) edgeSet.emplace(ed.first);
        if (!terminalsConnected(edgeSet, terminals)) {
            cerr << "Error: Unable to connect all terminals in the fallback solution.\n";
            exit(1);
        }
        best_init_cost = computeCost(edgeSet);
    }

    minVal = best_init_cost;
    bestSolutionEdges = best_init_solution;
    W = best_init_solution;

    // Apply enhanced local searches
    int max_iter = 50;
    for (int i = 0; i < max_iter; i++) {
        bool improved = local_search_improved(W, terminals, graph);
        if (improved) {
            // Update minVal and bestSolutionEdges are already handled in the function
            continue;
        }
        else {
            break; // No further improvements
        }
    }


    outputResults(W);

    return 0;
}
