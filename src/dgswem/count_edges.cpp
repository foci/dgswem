#include <vector>
#include <cstdio>

using std::vector;
using std::pair;

typedef int node;
typedef int edge;

vector<vector<pair<node, edge>>> Adj;

extern "C" {

  void create_adj_list(int np, int nedges, int* nedno1, int* nedno2) {
    Adj.resize(np);

    for (int i = 0; i < nedges; i++) {
      int n1 = nedno1[i] - 1;
      int n2 = nedno2[i] - 1;

      Adj[n1].push_back({n2, i+1});
      Adj[n2].push_back({n1, i+1});
    }

  }

  int get_edge_no(int n1, int n2) {
    int m1 = n1 - 1;
    int m2 = n2 - 1;
    for (int i = 0; i < Adj[m1].size(); i++) {
      if (Adj[m1][i].first == m2) {
        return Adj[m1][i].second;
      }
    }

    return 0;
  }




};
