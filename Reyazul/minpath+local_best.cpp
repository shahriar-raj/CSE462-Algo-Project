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

static long long minVal = LLONG_MAX;
static long long optimal_solution = 503; // Known optimal for demonstration
static vector<int> global_terminals;
static map<pair<int,int>, int> bestSolutionEdges;
static vector<bool> is_terminal;
static int n; // number of nodes

// The global edge map W stores all edges of the input graph in normalized form:
// For an edge (u,v), we always store it as (min(u,v), max(u,v)).
static map<pair<int,int>, int> original_W;

// Dijkstra's shortest path: returns prev[] for path reconstruction
vector<int> dijkstraShortestPath(const vector<vector<Edge>>& graph, int start) {
    int N = (int)graph.size();
    vector<int> dist(N, INT_MAX);
    vector<int> prev(N, -1);

    dist[start] = 0;
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
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

// Check if all terminals are connected in a given set of edges
bool terminalsConnected(const set<pair<int,int>>& edges, const vector<int>& terminals) {
    vector<vector<int>> adj(n+1);
    for (auto &ed : edges) {
        int u = ed.first, v = ed.second;
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    int start = terminals[0];
    unordered_set<int> visited;
    queue<int>q;
    q.push(start);
    visited.insert(start);

    while(!q.empty()) {
        int u = q.front(); q.pop();
        for (auto &nx : adj[u]) {
            if (!visited.count(nx)) {
                visited.insert(nx);
                q.push(nx);
            }
        }
    }

    for (auto term : terminals) {
        if (!visited.count(term)) return false;
    }
    return true;
}

// Compute cost of a solution given its edges and original_W
long long computeCost(const set<pair<int,int>>& edges) {
    long long cost = 0;
    for (auto &e : edges) {
        // e is already normalized
        cost += original_W.at(e);
    }
    return cost;
}

// Build a solution starting from a given terminal
map<pair<int,int>, int> buildSolutionFromTerminal(int startVertex,
                                                  const vector<int>& terminals,
                                                  const vector<vector<Edge>>& graph) {
    set<int> steinerTreeVertices;
    steinerTreeVertices.insert(startVertex);

    map<pair<int,int>, int> steinerTreeEdges;
    vector<int> remainingTerminals = terminals;
    remainingTerminals.erase(
        remove(remainingTerminals.begin(), remainingTerminals.end(), startVertex),
        remainingTerminals.end()
    );

    while (!remainingTerminals.empty()) {
        int nearestTerminal = -1;
        int minPathLength = INT_MAX;

        map<pair<int,int>, int> bestPathEdges;

        for (int treeVertex : steinerTreeVertices) {
            vector<int> prev = dijkstraShortestPath(graph, treeVertex);

            for (int terminal : remainingTerminals) {
                int current = terminal;
                int pathLength = 0;
                vector<pair<int,int>> pathEdges;

                while (current != -1 && current != treeVertex) {
                    int parent = prev[current];
                    if (parent == -1) break;

                    int a = min(parent, current);
                    int b = max(parent, current);
                    pair<int,int> ed(a,b);

                    // Ensure edge is in original_W
                    if (original_W.find(ed) == original_W.end()) {
                        // Edge not found, break this path
                        pathEdges.clear();
                        break;
                    }

                    pathLength += original_W.at(ed);
                    pathEdges.push_back(ed);
                    current = parent;
                }

                if (current == treeVertex && !pathEdges.empty() && pathLength < minPathLength) {
                    minPathLength = pathLength;
                    nearestTerminal = terminal;
                    bestPathEdges.clear();
                    for (auto edge : pathEdges) {
                        bestPathEdges[edge] = original_W.at(edge);
                    }
                }
            }
        }

        if (nearestTerminal == -1) {
            // No path found, break
            break;
        }

        for (auto& edge : bestPathEdges) {
            steinerTreeEdges[edge.first] = edge.second;
            steinerTreeVertices.insert(edge.first.first);
            steinerTreeVertices.insert(edge.first.second);
        }

        remainingTerminals.erase(
            remove(remainingTerminals.begin(), remainingTerminals.end(), nearestTerminal),
            remainingTerminals.end()
        );
    }

    return steinerTreeEdges;
}

// Local search: Remove non-terminal nodes if it reduces cost
bool local_search_remove_nodes(map<pair<int,int>, int>& W_map, const vector<int>& terminals) {
    set<pair<int,int>> edgeSet;
    for (auto &e : W_map) edgeSet.insert(e.first);

    long long best_cost = computeCost(edgeSet);

    // Identify nodes in solution
    unordered_set<int> nodes_in_solution;
    for (auto &ed : edgeSet) {
        nodes_in_solution.insert(ed.first);
        nodes_in_solution.insert(ed.second);
    }

    for (int node : nodes_in_solution) {
        if (is_terminal[node]) continue; // skip terminals

        // Attempt to remove this node
        set<pair<int,int>> new_solution = edgeSet;
        vector<pair<int,int>> to_remove;
        for (auto &ed : edgeSet) {
            if (ed.first == node || ed.second == node) {
                to_remove.push_back(ed);
            }
        }
        for (auto &r : to_remove) {
            new_solution.erase(r);
        }

        if (terminalsConnected(new_solution, terminals)) {
            long long new_cost = computeCost(new_solution);
            if (new_cost < best_cost) {
                // Found improvement
                map<pair<int,int>, int> newW;
                for (auto &ed : new_solution) {
                    newW[ed] = original_W.at(ed);
                }
                W_map = newW;
                minVal = new_cost;
                bestSolutionEdges = newW;
                return true;
            }
        }
    }

    return false;
}

// Local search: Try replacing edges
// For each edge in the solution, remove it and attempt to reconnect components cheaper.
bool local_search_edge_replacement(map<pair<int,int>, int>& W_map, const vector<int>& terminals, const vector<vector<Edge>>& graph) {
    set<pair<int,int>> edgeSet;
    for (auto &e : W_map) edgeSet.insert(e.first);
    long long best_cost = computeCost(edgeSet);

    // Try removing each edge and see if we can reconnect cheaper
    for (auto &removed_edge : edgeSet) {
        set<pair<int,int>> new_solution = edgeSet;
        new_solution.erase(removed_edge);

        if (terminalsConnected(new_solution, terminals)) {
            // Cheaper by just removing edge?
            long long new_cost = computeCost(new_solution);
            if (new_cost < best_cost) {
                map<pair<int,int>, int> newW;
                for (auto &ed : new_solution) {
                    newW[ed] = original_W.at(ed);
                }
                W_map = newW;
                minVal = new_cost;
                bestSolutionEdges = newW;
                return true;
            }
            continue;
        }

        // Not connected after removal, we must reconnect:
        // Identify two components:
        vector<vector<int>> adj(n+1);
        for (auto &ed : new_solution) {
            adj[ed.first].push_back(ed.second);
            adj[ed.second].push_back(ed.first);
        }

        // BFS from terminals[0] to find component A
        unordered_set<int> compA;
        {
            queue<int>q;
            q.push(terminals[0]);
            compA.insert(terminals[0]);
            while(!q.empty()) {
                int u = q.front(); q.pop();
                for (auto &nx : adj[u]) {
                    if (!compA.count(nx)) {
                        compA.insert(nx);
                        q.push(nx);
                    }
                }
            }
        }

        // Component B = nodes_in_solution - compA
        unordered_set<int> nodes_in_solution_set;
        for (auto &ed : new_solution) {
            nodes_in_solution_set.insert(ed.first);
            nodes_in_solution_set.insert(ed.second);
        }

        unordered_set<int> compB;
        for (auto node : nodes_in_solution_set) {
            if (!compA.count(node)) compB.insert(node);
        }

        if (compB.empty()) continue;

        // Try to find a cheaper path connecting compA to compB
        // We'll run Dijkstra from each node in compA and look for shortest path to compB
        long long best_path_cost = LLONG_MAX;
        vector<pair<int,int>> best_path_edges_to_add;

        for (int source : compA) {
            vector<int> dist(n+1, INT_MAX), parent(n+1, -1);
            dist[source] = 0;
            priority_queue<pair<int,int>, vector<pair<int,int>>, greater<>> pq;
            pq.push({0, source});

            while(!pq.empty()) {
                auto [cd, u] = pq.top(); pq.pop();
                if (cd > dist[u]) continue;
                for (auto &ed : graph[u]) {
                    int v = ed.to;
                    int nd = cd + ed.length;
                    if (nd < dist[v]) {
                        dist[v] = nd;
                        parent[v] = u;
                        pq.push({nd,v});
                    }
                }
            }

            for (int target : compB) {
                if (dist[target] == INT_MAX) continue;
                // Reconstruct path
                vector<int> path;
                int cur = target;
                while (cur != -1) {
                    path.push_back(cur);
                    cur = parent[cur];
                }
                reverse(path.begin(), path.end());

                // Convert path to edges and compute cost
                long long path_cost = 0;
                vector<pair<int,int>> path_edges;
                bool valid_path = true;
                for (size_t i = 0; i+1 < path.size(); i++) {
                    int a = min(path[i], path[i+1]);
                    int b = max(path[i], path[i+1]);
                    pair<int,int> ped(a,b);
                    if (original_W.find(ped) == original_W.end()) {
                        // edge not found, invalid path
                        valid_path = false;
                        break;
                    }
                    path_cost += original_W.at(ped);
                    path_edges.push_back(ped);
                }

                if (valid_path && path_cost < best_path_cost) {
                    best_path_cost = path_cost;
                    best_path_edges_to_add = path_edges;
                }
            }
        }

        if (best_path_cost < LLONG_MAX) {
            // new cost = cost of new_solution without removed_edge + best_path_cost
            long long new_cost_without = computeCost(new_solution);
            long long new_cost = new_cost_without + best_path_cost;
            if (new_cost < best_cost) {
                // Update solution
                map<pair<int,int>, int> newWmap;
                for (auto &ed : new_solution) {
                    newWmap[ed] = original_W.at(ed);
                }
                for (auto &pe : best_path_edges_to_add) {
                    newWmap[pe] = original_W.at(pe);
                }
                // Check connectivity again for safety
                set<pair<int,int>> finalEdges;
                for (auto &x : newWmap) finalEdges.insert(x.first);
                if (terminalsConnected(finalEdges, terminals)) {
                    W_map = newWmap;
                    minVal = new_cost;
                    bestSolutionEdges = newWmap;
                    return true;
                }
            }
        }
    }

    return false;
}

void parseGraphSection(vector<vector<Edge>>& graph, map<pair<int,int>, int>& W,
                       vector<long>& incident, set<long>& weightSet) {
    string dummy;
    long m, u, v, w;
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
        pair<int,int> ed(min(u,v), max(u,v));
        W[ed] = w;
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
        // Process if needed
    }

    long tu, tv;
    for (long i = 0; i < b - 1; i++) {
        cin >> tu >> tv;
    }
    cin >> dummy; // end
}

void outputResults(const map<pair<int,int>, int>& W_map) {
    cout << "VALUE " << minVal << endl;
    if (optimal_solution > 0) {
        double approx = (double)minVal / (double)optimal_solution;
        cout << "APPROXIMATION RATIO " << approx << "\n";
    }

    cout << "EDGES" << endl;
    for (const auto& edge : W_map) {
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

    original_W = W; // Keep a global reference to all edges
    global_terminals = terminals;
    is_terminal.assign(n+1, false);
    for (auto &term : terminals) {
        is_terminal[term] = true;
    }

    // Try multiple initial solutions from different terminals
    long long best_init_cost = LLONG_MAX;
    map<pair<int,int>, int> best_init_solution;

    for (auto startT : terminals) {
        auto candidate = buildSolutionFromTerminal(startT, terminals, graph);
        set<pair<int,int>> candEdges;
        for (auto &ed : candidate) candEdges.insert(ed.first);
        if (terminalsConnected(candEdges, terminals)) {
            long long cost = computeCost(candEdges);
            if (cost < best_init_cost) {
                best_init_cost = cost;
                best_init_solution = candidate;
            }
        }
    }

    if (best_init_solution.empty()) {
        // fallback if no solution found
        best_init_solution = buildSolutionFromTerminal(terminals[0], terminals, graph);
        set<pair<int,int>> edgeSet;
        for (auto &ed : best_init_solution) edgeSet.insert(ed.first);
        best_init_cost = computeCost(edgeSet);
    }

    minVal = best_init_cost;
    bestSolutionEdges = best_init_solution;
    W = best_init_solution;

    // Apply local searches
    int max_iter = 50;
    for (int i = 0; i < max_iter; i++) {
        bool improved = false;
        improved = local_search_remove_nodes(W, terminals);
        if (improved) continue;
        improved = local_search_edge_replacement(W, terminals, graph);
        if (!improved) break;
    }

    outputResults(W);

    return 0;
}
