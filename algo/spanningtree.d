module algo.spanningtree;

private import jive.array;
private import jive.priorityqueue;
private import jive.unionfind;

import algo.graph;

/**
 * compute a minimum spanning tree (or forest if the graph is disconnected)
 * implemented using Kruskal's algorithm
 */
Graph!(false, T) minimumSpanningTree(G, T = G.Label)(G g)
{
    auto tree = Graph!(false, T)(g.n);

    Array!(Edge!T) edges;
    for(int a = 0; a < g.n; ++a)
        foreach(b, c; g.succ(a))
            edges.pushBack(Edge!T(a, b, c));
    sort(edges[]);

    auto uf = UnionFind(g.n);
    foreach(e; edges[])
        if(!uf.isJoined(e.from, e.to))
        {
            tree.addEdge(e.from, e.to, e.label);
            uf.join(e.from, e.to);
        }

    return tree;
}
