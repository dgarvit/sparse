#ifndef MAKE_SET_H
#define MAKE_SET_H

#include <set>
#include <vector>

void visit(std::vector<std::set<int>> const &dag, int idx, bool visited[], std::set<int>& reachset) {
    reachset.insert(idx);
    visited[idx] = true;
    for (auto f : dag[idx]) {
        if (!visited[f])
            visit(dag, f, visited, reachset);
    }
}

void insertToLevel(std::vector< std::vector<int> >& levels, int const node, int const currLevel) {
    if (levels.size() <= currLevel) {
        levels.emplace_back(); // construct a new level
    }
    levels[currLevel].push_back(node);
}

// Function that returns the DAG. visited matrix is later used to generate reachset
std::vector<std::set<int>> makeDAG(int* &col, int* &row, double* &val, int n, bool* &visited) {
    std::vector<std::set<int>> dag(n);
    for (int i = 0; i < n; ++i) {
        visited[i] = false;
    }
    

    // Making the dag
    for (int i = 1; i < n; ++i) {
        for (int j = col[i-1]; j < col[i]; ++j) {
            if (row[j] != i - 1) // Do not create a self-loop
                dag[i-1].insert(row[j]);
        }
    }
    return dag;
}

// Function that generates the reahset from the DAG
std::set<int> makeReachset(std::vector<std::set<int>> const &dag, bool* &visited, int n, double* &x) {
    std::set<int> reachset;

    // Creating dependencies using DFS
    for (int i = 0; i < n; ++i) {
        if (x[i] != 0) {
            if (!visited[i]) {
                visit(dag, i, visited, reachset);
            }
        }
    }
    return reachset;
}

// Creating levels and tasks
std::vector< std::vector<int> > makeLevelSet(std::vector<std::set<int>> const &dag, std::set<int> reachset, int n) {
    std::vector<int> indegree(n); // Keep a count of incoming edges to record root nodes
    std::vector< std::vector<int> > levels; // Store levels of parallelism. Tasks inside same level can be performed parallely; each 2 levels need to sync
    int currLevel = 0;

    for (auto f : reachset) {
        for (auto j : dag[f]) {
            ++indegree[j];
        }
    }

    while (!reachset.empty()) {
        for (auto f : reachset) {
            if (indegree[f] == 0) {
                insertToLevel(levels, f, currLevel);
            }
        }

        for (auto f : levels[currLevel]) {
            reachset.erase(f);
            for (auto j : dag[f]) {
                --indegree[j];
            }
        }
        ++currLevel;
    }
    return levels;
}

#endif
