#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <algorithm>
#include <climits>
#include <cmath>
#include <set>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

class SteinerTree {
private:
    int n, e, t;  // nodes, edges, terminals
    long long max_cost;
    long long optimal_solution;  // Hardcoded optimal solution

    // Graph representation: adjacency list {node -> {neighbor -> weight}}
    unordered_map<int, unordered_map<int, int>> graph;
    unordered_map<int, bool> is_terminal;

    vector<int> terminals;

    // For Dijkstra
    struct DNode {
        int node;
        long long dist;
        bool operator>(const DNode &other) const {
            return dist > other.dist;
        }
    };

    // For MST
    struct Edge {
        int u, v;
        long long w;
    };

    vector<vector<long long>> all_dist;
    vector<vector<int>> all_parent; // all_parent[i][v] = parent of v in shortest path tree from terminals[i]

    // Steiner tree solution representation: set of edges (u < v)
    set<pair<int,int>> steiner_edges;

    void read_input(const string& filename) {
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "Unable to open file: " << filename << endl;
            exit(1);
        }

        string line, section;
        max_cost = 0;

        while (getline(file, line)) {
            istringstream iss(line);
            string token;
            iss >> token;

            if (token == "SECTION") {
                iss >> section;

                if (section == "Graph") {
                    // Parse graph section
                    while (getline(file, line) && line.find("END") == string::npos) {
                        iss.clear();
                        iss.str(line);
                        string type;
                        iss >> type;

                        if (type == "Nodes") {
                            iss >> n;
                        } else if (type == "Edges") {
                            iss >> e;
                        } else if (type == "E") {
                            int u, v, w;
                            iss >> u >> v >> w;

                            graph[u][v] = w;
                            graph[v][u] = w;
                            max_cost += w;
                            if (!is_terminal.count(u)) is_terminal[u] = false;
                            if (!is_terminal.count(v)) is_terminal[v] = false;
                        }
                    }
                }
                else if (section == "Terminals") {
                    // Parse terminals section
                    while (getline(file, line) && line.find("END") == string::npos) {
                        iss.clear();
                        iss.str(line);
                        string type;
                        iss >> type;

                        if (type == "Terminals") {
                            iss >> t;
                        } else if (type == "T") {
                            int terminal;
                            iss >> terminal;
                            terminals.push_back(terminal);
                            is_terminal[terminal] = true;
                        }
                    }
                }
            }

            if (line.find("EOF") != string::npos) {
                break;
            }
        }

        file.close();

        // Sort terminals
        sort(terminals.begin(), terminals.end());
    }

    void dijkstra(int start, vector<long long>& dist, vector<int>& parent) {
        dist.assign(n+1, LLONG_MAX);
        parent.assign(n+1, -1);
        dist[start] = 0;
        priority_queue<DNode, vector<DNode>, greater<DNode>> pq;
        pq.push({start, 0});

        while(!pq.empty()) {
            auto [node, d] = pq.top();
            pq.pop();
            if (d > dist[node]) continue;
            for (auto &ed : graph[node]) {
                int nxt = ed.first; int w = ed.second;
                if (dist[node] + w < dist[nxt]) {
                    dist[nxt] = dist[node] + w;
                    parent[nxt] = node;
                    pq.push({nxt, dist[nxt]});
                }
            }
        }
    }

    void compute_terminal_distances() {
        int T = (int)terminals.size();
        all_dist.assign(T, vector<long long>(n+1, LLONG_MAX));
        all_parent.assign(T, vector<int>(n+1, -1));

        for (int i = 0; i < T; i++) {
            int start = terminals[i];
            dijkstra(start, all_dist[i], all_parent[i]);
        }
    }

    long long get_term_dist(int i, int j) {
        // shortest distance between terminals[i] and terminals[j]
        return all_dist[i][terminals[j]];
    }

    vector<int> reconstruct_path(int i, int j) {
        // Reconstruct path between terminals[i] and terminals[j]
        int start = terminals[i];
        int goal = terminals[j];
        vector<int> path;
        int cur = goal;
        while (cur != -1 && cur != start) {
            path.push_back(cur);
            cur = all_parent[i][cur];
        }
        path.push_back(start);
        reverse(path.begin(), path.end());
        return path;
    }

    vector<pair<int,int>> compute_metric_mst() {
        int T = (int)terminals.size();
        // Build metric closure MST using Prim's
        vector<bool> in_mst(T, false);
        vector<long long> key(T, LLONG_MAX);
        vector<int> parent(T, -1);

        key[0] = 0;
        for (int count = 0; count < T-1; count++) {
            int u = -1; long long mn = LLONG_MAX;
            for (int v = 0; v < T; v++) {
                if (!in_mst[v] && key[v] < mn) {
                    mn = key[v];
                    u = v;
                }
            }
            in_mst[u] = true;
            for (int w = 0; w < T; w++) {
                if (!in_mst[w]) {
                    long long d = get_term_dist(u, w);
                    if (d < key[w]) {
                        key[w] = d;
                        parent[w] = u;
                    }
                }
            }
        }

        vector<pair<int,int>> mst_edges;
        for (int i = 1; i < T; i++) {
            mst_edges.push_back({parent[i], i});
        }
        return mst_edges;
    }

    void build_steiner_tree(const vector<pair<int,int>>& mst_edges) {
        steiner_edges.clear();
        for (auto &ed : mst_edges) {
            int i = ed.first; int j = ed.second;
            vector<int> path = reconstruct_path(i, j);
            for (size_t k = 0; k+1 < path.size(); k++) {
                int u = path[k]; int v = path[k+1];
                if (u > v) std::swap(u, v);
                steiner_edges.insert({u, v});
            }
        }
    }

    long long compute_steiner_cost(const set<pair<int,int>>& edges) {
        long long total = 0;
        for (auto &ed : edges) {
            int u = ed.first, v = ed.second;
            total += graph[u][v];
        }
        return total;
    }

    long long current_cost() {
        return compute_steiner_cost(steiner_edges);
    }

    bool terminals_connected(const set<pair<int,int>>& edges) {
        // Check if all terminals are in one connected component
        // We'll build adjacency from edges and run BFS or Union-Find
        unordered_map<int, vector<int>> adj;
        for (auto &ed : edges) {
            int u = ed.first; int v = ed.second;
            adj[u].push_back(v);
            adj[v].push_back(u);
        }

        // BFS from first terminal
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

        // Check all terminals in visited
        for (auto term : terminals) {
            if (!visited.count(term)) return false;
        }

        return true;
    }

    // Local Search:
    // Try removing non-terminal nodes if possible and beneficial
    // A non-terminal node is any node in the steiner solution that is not a terminal.
    // If removing it (and its incident edges) still keeps all terminals connected and lowers the cost, accept it.
    bool local_search_improve() {
        // Build adjacency from current steiner_edges
        // Identify nodes in the solution
        unordered_set<int> nodes_in_solution;
        for (auto &ed : steiner_edges) {
            nodes_in_solution.insert(ed.first);
            nodes_in_solution.insert(ed.second);
        }

        long long best_cost = current_cost();
        set<pair<int,int>> best_edges = steiner_edges;

        // Try removing each non-terminal node that is in the solution
        for (int node : nodes_in_solution) {
            if (is_terminal[node]) continue; // skip terminals

            // Remove all edges incident to 'node'
            set<pair<int,int>> new_solution = steiner_edges;
            // Collect edges to remove
            vector<pair<int,int>> to_remove;
            for (auto &ed : steiner_edges) {
                int u = ed.first; int v = ed.second;
                if (u == node || v == node) {
                    to_remove.push_back(ed);
                }
            }

            for (auto &r : to_remove) {
                new_solution.erase(r);
            }

            // Check feasibility
            if (terminals_connected(new_solution)) {
                long long new_cost = compute_steiner_cost(new_solution);
                if (new_cost < best_cost) {
                    // Found improvement
                    steiner_edges = new_solution;
                    return true;
                }
            }
        }

        // If no improvement found
        return false;
    }

public:
    SteinerTree(long long optimal_sol = 0) :
        max_cost(0), n(0), e(0), t(0), optimal_solution(optimal_sol) {}

    void solve_from_file(const string& filename) {
        read_input(filename);
        solve();
    }

    void solve() {
        // Step 1: Compute shortest paths from each terminal
        compute_terminal_distances();

        // Step 2: Build metric closure MST on terminals
        auto mst_edges = compute_metric_mst();

        // Step 3: Reconstruct Steiner solution from MST
        build_steiner_tree(mst_edges);

        // Step 4: Apply local search improvements
        // We'll do a simple loop until no improvement or max iterations
        int max_iter = 50; // arbitrary small limit
        for (int i = 0; i < max_iter; i++) {
            bool improved = local_search_improve();
            if (!improved) break;
        }

        long long best_solution = current_cost();

        cout << "VALUE " << best_solution << "\n";
        if (optimal_solution > 0) {
            double approx_ratio = (double)best_solution / (double)optimal_solution;
            cout << "APPROXIMATION RATIO " << approx_ratio << "\n";
        }
    }
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);

    // Hardcoded input file name and optimal solution for demonstration
    string input_file = "./instance025.gr";
    long long optimal_solution = 5194; // Replace with actual optimal solution if known

    SteinerTree steiner(optimal_solution);
    steiner.solve_from_file(input_file);

    return 0;
}
