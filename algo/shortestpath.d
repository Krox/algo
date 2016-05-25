module algo.shortestpath;

/**

algorithms belonging in this module

    single-source without negative edges:
        Dijkstra with naive selection   O(V^2)                  -
        Dijkstra with binary heap       O(E * log(V))           done
        Dijkstra with Fibonacci heap    O(V * log(V) + E)       -

    single-source with negative edges:
        Bellman–Ford                    O(V * E)                TODO

    all-pairs with/without negative edges:
        Floyd-Warshall                  O(V^3)                  done
        Johnson with binary heap        O(V * E * log(V))       TODO
        Johnson with Fibonacci heap     O(V * E + V^2 * log(V)) -

**/

private import std.algorithm;
private import jive.array;
private import jive.bitarray;
private import jive.priorityqueue;

import algo.graph;

/**
 * compute shortest path from one vertex to every other vertex. O(E * log(V)).
 */
struct Dijkstra(T)
{
    // start vertex has pred = -1 and dist = 0
    // unreachable vertices have pred = -1 and dist = T.max

    Array!T dist;
    Array!(HalfEdge!T) pred;

    @disable this();
    @disable this(this);

    this(G)(auto ref G g, int start)
    {
        dist.resize(g.n, T.max);
        pred.resize(g.n, HalfEdge!T(-1, T.init));
        PriorityQueue!(HalfEdge!T) q;

        q.push(HalfEdge!T(start, 0));
        dist[start] = 0;

        while(!q.empty)
        {
            auto e = q.pop;

            if(e.label > dist[e.to])
                continue;
            assert(e.label == dist[e.to]);

            foreach(b, c; g.succ(e.to))
                if(e.label + c < dist[b])
                {
                    dist[b] = e.label + c;
                    q.push(HalfEdge!T(b, e.label + c));
                    pred[b] = HalfEdge!T(e.to, c);
                }
        }
    }

    /** list all edges of the shortest-path-tree */
    Edge!T[] listEdges()
    {
        Array!(Edge!T) r;
        foreach(int x, HalfEdge!T e; pred)
            if(e.to != -1)
                r.pushBack(Edge!T(e.to, x, e.label));
        return r.release();
    }

    /** list all the edges of the shortest path leading to vertex v */
    Edge!T[] listEdges(int v)
    {
        Array!(Edge!T) r;
        while(pred[v].to != -1)
        {
            r.pushBack(Edge!T(pred[v].to, v, pred[v].label));
            v = pred[v].to;
        }
        reverse(r[]);
        return r.release();
    }
}

/** convenience wrapper for better template parameter deduction */
Dijkstra!T dijkstra(G, T = G.Label)(G g, int start)
{
    return Dijkstra!T(g, start);
}

/**
 * compute shortest path from every vertex to every other vertex. O(V^3).
 */
struct FloydWarshall(T)
{
	Array2!T dist;
    Array2!(HalfEdge!T) pred;

	@disable this();
	@disable this(this);

	this(G)(auto ref G g)
	{
		dist = Array2!T(g.n, g.n, T.max);
        pred = Array2!(HalfEdge!T)(g.n, g.n, HalfEdge!T(-1, T.init));

        for(int i = 0; i < g.n; ++i)
            dist[i,i] = 0;

		for(int i = 0; i < g.n; ++i)
			foreach(j, c; g.succ(i))
            {
				dist[i,j] = c;
                pred[i,j] = HalfEdge!T(i, c);
            }

		for(int b = 0; b < g.n; ++b)
			for(int a = 0; a < g.n; ++a)
				if(dist[a,b] < T.max)
					for(int c = 0; c < g.n; ++c)
                        if(dist[b,c] < T.max && dist[a,b] + dist[b,c] < dist[a,c])
                        {
					        dist[a,c] = dist[a,b] + dist[b,c];
                            pred[a,c] = pred[b,c];
                        }
	}

    /** list all edges of the shortest-path-tree starting at vertex v */
    Edge!T[] listEdges(int v)
    {
        Array!(Edge!T) r;
        for(int w = 0; w < dist[].size[0]; ++w)
            if(pred[v,w].to != -1)
                r.pushBack(Edge!T(pred[v,w].to, w, pred[v,w].label));
        return r.release();
    }
}

/** convenience wrapper for better template parameter deduction */
FloydWarshall!T floydWarshall(G, T = G.Label)(G g)
{
    return FloydWarshall!T(g);
}
