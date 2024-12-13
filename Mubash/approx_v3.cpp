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
    
    // Graph representation
    unordered_map<int, unordered_map<int, int>> graph;
    unordered_map<int, bool> is_terminal;
    
    // Steiner tree components
    unordered_map<int, unordered_map<int, int>> steiner_tree;
    unordered_map<int, unordered_map<int, int>> steiner_cost;
    
    vector<int> terminals;

    // Helper structs for priority queue
    struct QueueNode {
        int node, distance;
        vector<int> path;
        bool operator>(const QueueNode& other) const {
            return distance > other.distance;
        }
    };

    // Helpers for set operations and combinations
    vector<vector<int>> get_all_subsets(const vector<int>& source) {
        vector<vector<int>> subsets;
        int n = source.size();
        
        // Generate all possible subset combinations
        for (int i = 1; i < (1 << n); ++i) {
            vector<int> subset;
            for (int j = 0; j < n; ++j) {
                if (i & (1 << j)) {
                    subset.push_back(source[j]);
                }
            }
            subsets.push_back(subset);
        }
        return subsets;
    }

    // Dijkstra-based shortest path 
    unordered_map<int, unordered_map<int, int>> dijkstra(int start) {
        unordered_map<int, unordered_map<int, int>> distances;
        priority_queue<QueueNode, vector<QueueNode>, greater<QueueNode>> pq;
        
        pq.push({start, 0, {start}});
        
        while (!pq.empty()) {
            QueueNode current = pq.top();
            pq.pop();
            
            int node = current.node;
            int dist = current.distance;
            
            // Skip if already processed
            if (distances[start].count(node) && distances[start][node] <= dist) 
                continue;
            
            distances[start][node] = dist;
            
            // Explore neighbors
            for (auto& [neighbor, edge_cost] : graph[node]) {
                int new_dist = dist + edge_cost;
                
                if (!distances[start].count(neighbor) || new_dist < distances[start][neighbor]) {
                    pq.push({neighbor, new_dist, current.path});
                }
            }
        }
        
        return distances;
    }

    // Tree cost calculation
    long long calculate_tree_cost() {
        long long total_cost = 0;
        set<pair<int, int>> unique_edges;
        
        for (auto& [u, neighbors] : steiner_tree) {
            for (auto& [v, cost] : neighbors) {
                unique_edges.insert({min(u, v), max(u, v)});
            }
        }
        
        for (auto& [u, v] : unique_edges) {
            total_cost += graph[u][v];
        }
        
        return total_cost;
    }

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
                            is_terminal[u] = false;
                            is_terminal[v] = false;
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

public:
    SteinerTree(long long optimal_sol = 0) : 
        max_cost(0), n(0), e(0), t(0), optimal_solution(optimal_sol) {}

    void solve_from_file(const string& filename) {
        read_input(filename);
        solve();
    }

    void solve() {
        long long best_solution = max_cost + 1;
        set<pair<int, int>> best_tree;
        
        int limit = min(t, 200);
        vector<int> check_terminals(terminals.begin(), terminals.begin() + limit);
        
        while (!check_terminals.empty()) {
            int start_terminal = check_terminals.front();
            check_terminals.erase(check_terminals.begin());
            
            // Reset Steiner tree
            steiner_tree.clear();
            steiner_cost = graph;
            
            // Compute distances from start terminal
            auto distances = dijkstra(start_terminal);
            
            // Sort other terminals by distance
            vector<int> other_terminals = terminals;
            other_terminals.erase(
                remove(other_terminals.begin(), other_terminals.end(), start_terminal), 
                other_terminals.end()
            );
            
            sort(other_terminals.begin(), other_terminals.end(), 
                [&](int a, int b) { return distances[start_terminal].at(a) < distances[start_terminal].at(b); });
            
            // Connect first terminal
            int closest_terminal = other_terminals[0];
            
            // Rest of the optimization logic would go here...
        }
        
        // Output result with approximation ratio
        cout << "VALUE " << best_solution << endl;
        
        // Calculate and print approximation ratio if optimal solution is known
        if (optimal_solution > 0) {
            double approx_ratio = static_cast<double>(best_solution) / optimal_solution;
            cout << "APPROXIMATION RATIO " << fixed << approx_ratio << endl;
        }
    }
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);

    // Hardcoded input file name and optimal solution
    string input_file = "./Track3/instance001.gr";
    long long optimal_solution = 2256; // Replace with actual optimal solution

    SteinerTree steiner(optimal_solution);
    steiner.solve_from_file(input_file);

    return 0;
}