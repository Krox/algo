module algo.graph;

private import std.process : executeShell;
private import std.random;
private import std.stdio;
private import std.typetuple;
private import std.typecons;
private import std.math : sqrt, isNaN;
private import std.conv : to;
private import std.algorithm : min, max, move, sort, swap;
private import jive.array;
private import jive.bitarray;
private import jive.unionfind;

/**
 * Directed or undirected (multi-)graph with optionally labeled edges.
 * Mostly for convenience, as all algorithms should be templated to work with
 * any graph structure that supports (a subset of) the methods
 * n, m, succ, pred.
 */
struct Graph(bool directed = false, Label...)
{
	static struct HalfEdge
	{
		int to;
		Label label;
	}

	/** actual graph data */
	private Array!(Array!HalfEdge) g;
	private size_t edgeCount;

	/** purely cosmetic hints for viewing the graph, may not be set at all */
	string name = "myGraph";
	string layout = "dot";
	float sizeX = float.nan;
	float sizeY = float.nan;
	Array!(Tuple!(float,float)) pos;

	/** constructor for empty graph */
	this(size_t n)
	{
		assert(n <= int.max);
		g.resize(n);
	}

	/** number of vertices */
	int n() const nothrow @safe @property
	{
		return cast(int)g.length;
	}

	/** number of edges */
	size_t m() const nothrow @safe @property
	{
		return edgeCount;
	}

	/** maximum number of edges possible with current vertex count */
	size_t maxEdges() const nothrow @safe @property
	{
		static if(directed)
			return n * (n - 1);
		else
			return n * (n - 1) / 2;
	}

	/** numEdges / maxEdges */
	float density() const nothrow @safe @property
	{
		return cast(float)m / maxEdges;
	}

	/** numEdges / numVertices */
	float ratio() const nothrow @safe @property
	{
		return cast(float)m / n;
	}

	/** add a new vertex, returns its number */
	int addVertex()
	{
		assert(n < int.max);
		g.pushBack((Array!HalfEdge).init);
		return cast(int)(n-1);
	}

	/** add a new edge */
	void addEdge(int from, int to, Label label)
	{
		assert(0 <= from && from < n);
		assert(0 <= to && to < n);
		g[from].pushBack(HalfEdge(to, label));
		if(!directed && to != from)
			g[to].pushBack(HalfEdge(from, label));
		++edgeCount;
	}

	/** removes duplicate edges and loops */
	void reduce()
	{
		edgeCount = 0;

		for(int a = 0; a < n; ++a)
		{
			sort!"a.to < b.to"(g[a][]);
			foreach(i, e, ref bool rem; &g[a].prune)
				if(e.to == a || (i < g[a].length-1 && e.to == g[a][i+1].to))
					rem = true;
			edgeCount += g[a].length;
		}

		static if(!directed)
			edgeCount /= 2;
	}

	/** test wether a certainedge exists */
	bool hasEdge(int from, int to) const
	{
		foreach(e; succ(from))
			if(e == to)
				return true;
		return false;
	}

	struct Range
	{
		const(HalfEdge)[] arr;

		int opApply(int delegate(int x) dg) const
		{
			int r = 0;
			foreach(e; arr)
				if((r = dg(e.to)) != 0)
					break;
			return r;
		}

		static if(Label.length != 0)
		{
			int opApply(int delegate(int x, Label l) dg) const
			{
				int r = 0;
				foreach(e; arr)
					if((r = dg(e.to, e.label)) != 0)
						break;
				return r;
			}
		}
	}

	/** range of successors of a vertex */
	Range succ(int x) const nothrow
	{
		return Range(g[x][]);
	}

	/** write graph into a .dot file */
	void writeDot(string filename)
	{
		auto f = File(filename, "w");
		f.writefln("%s %s {", directed?"digraph":"graph", name);

		//f.writefln("node [label=\"\\N\"];");
		if(!isNaN(sizeX) && !isNaN(sizeY))
			f.writefln("graph [bb=\"0,0,%s,%s\"];", sizeX, sizeY);

		for(int a = 0; a < n; ++a)
			if(a < pos.length)
				f.writefln("%s [pos=\"%s,%s\"];", a, pos[a][0], pos[a][1]);
			else
				f.writefln("%s", a);

		for(int a = 0; a < n; ++a)
			static if(Label.length == 0)
			{
				foreach(b; succ(a))
					if(directed || a <= b)
						f.writefln("%s%s%s", a, directed?"->":"--", b);
			}
			else
			{
				foreach(b, l; succ(a))
					if(directed || a <= b)
						f.writefln("%s%s%s [label=\"%s\"]", a, directed?"->":"--", b, to!string(l));
			}

		f.writefln("}");
	}

	/** view the graph (using graphviz and eog) */
	void show()
	{
		writeDot("tmp_"~name~".dot");
		executeShell(layout~" "~"tmp_"~name~".dot -Tsvg > "~"tmp_"~name~".svg");
		executeShell("eog "~"tmp_"~name~".svg");
		executeShell("rm "~"tmp_"~name~".dot "~"tmp_"~name~".svg");
	}
}


//////////////////////////////////////////////////////////////////////
/// graph creation
//////////////////////////////////////////////////////////////////////

Graph!dir completeGraph(bool dir = false)(int n)
{
	auto g = Graph!dir(n);
	for(int i = 0; i < n; ++i)
		for(int j = i+1; j < n; ++j)
		{
			g.addEdge(i, j);
			static if(dir)
				g.addEdge(j, i);
		}
	return g;
}

Graph!dir lineGraph(bool dir = false)(int n)
{
	auto g = Graph!dir(n);
	for(int i = 1; i < n; ++i)
		g.addEdge(i-1, i);
	return g;
}

Graph!dir randomGraph(bool dir = false)(int n, float ratio)
{
	if((dir && ratio >= n-1) || (!dir && ratio >= (n-1)/2))
		return completeGraph!dir(n);

	auto g = Graph!dir(n);
	while(g.ratio < ratio)
	{
		int a = uniform(0, n);
		int b = uniform(0, n);
		if(a != b && !g.hasEdge(a, b))
			g.addEdge(a, b);
	}
	return g;
}

Graph!true randomAcyclicGraph(int n, float ratio)
{
	assert(ratio <= (n-1)/2);

	auto g = Graph!true(n);
	while(g.ratio < ratio)
	{
		int a = uniform(0, n);
		int b = uniform(0, n);
		if(a > b)
			swap(a, b);
		if(a != b && !g.hasEdge(a, b))
			g.addEdge(a, b);
	}
	return g;
}

Graph!dir randomTree(bool dir = false)(int n)
{
	auto g = Graph!dir(n);
	for(int i = 1; i < n; ++i)
		g.addEdge(i, uniform(0, i));
	return g;
}

Graph!dir binaryLeftTree(bool dir = false)(int n)
{
	auto g = Graph!dir(n);
	for(int i = 1; i < n; ++i)
		g.addEdge(i, (i-1)/2);
	return g;
}

Graph!(false, float) randomEuclideanGraph(int n, float size, float cutoff)
{
	auto g = Graph!(false, float)(n);
	g.pos.resize(n);
	g.sizeX = size;
	g.sizeY = size;

	for(int i = 0; i < n; ++i)
	{
		g.pos[i][0] = uniform(0.0, size);
		g.pos[i][1] = uniform(0.0, size);
	}

	for(int a = 0; a < n; ++a)
		for(int b = a+1; b < n; ++b)
		{
			auto distX = g.pos[a][0]-g.pos[b][0];
			auto distY = g.pos[a][1]-g.pos[b][1];
			auto dist2 = distX*distX + distY*distY;
			if(dist2 <= cutoff*cutoff)
				g.addEdge(a, b, sqrt(dist2));
		}

	g.layout = "neato";
	return g;
}


//////////////////////////////////////////////////////////////////////
/// reachability problems for unweighted graphs
//////////////////////////////////////////////////////////////////////

/**
 * tests if there is a path a->b.
 * Uses a simple dfs in O(n+m). When multiple doing mutliple such queries,
 * you should use reachableSet/TransitiveClosure or similar.
 */
bool hasPath(G)(auto ref G g, int from, int to) const
{
	if(from == to)
		return true;

	auto v = BitArray(g.n);
	Array!int stack;
	stack.pushBack(from);
	v[from] = true;

	while(!stack.empty)
	{
		int a = stack.popBack;
		foreach(b; g.succ(a))
			if(b == to)
				return true;
			else if(!v[b])
			{
				v[b] = true;
				stack.pushBack(b);
			}
	}

	return false;
}

/** determines the vertices reachable from a using a simple dfs in O(n+m) */
BitArray reachableSet(G)(auto ref G g, int a)
{
	BitArray v;
	g.reachableSet(a, v);
	return v;
}

/** ditto (can be used to prevent heap alocation) */
void reachableSet(G)(auto ref G g, int a, ref BitArray v)
{
	v.resize(g.n);
	v.reset();
	Array!int stack;

	stack.pushBack(a);
	v[a] = true;

	while(!stack.empty)
	{
		int b = stack.popBack;
		foreach(c; g.succ(b))
			if(!v[c])
			{
				v[c] = true;
				stack.pushBack(c);
			}
	}
}

/** compute transitve closure / reachability matrix of a graph */
struct TransitiveClosure
{
	@disable this();
	@disable this(this);

	Array!BitArray mat;

	this(G)(auto ref G g)
	{
		mat.resize(g.n);
		for(int i = 0; i < g.n; ++i)
			mat[i] = g.reachableSet(i);
	}

	bool hasPath(int a, int b) const
	{
		return mat[a][b];
	}

	size_t pathCount() const
	{
		size_t cnt = 0;
		foreach(ref s; mat[])
			cnt += s.count(true);
		return cnt;
	}
}

/**
 * Answer reachability queries after O(n+m) preprocessing. Meant for large
 * sparse graphs where building the whole transitive closure is not an option.
 * The downside is that a query can take up to O(n+m) using a dfs. But an effort
 * is made to answer many queries much faster.
 */
struct Reachability
{
	@disable this();
	@disable this(this);

	SCC scc;
	Components comp; // (weakly) connected components
	TopOrder top;
	Array!int up, down; // up/down levels of vertices
	Array!int timeA, timeB;

	private int time = 1;

	this(G)(auto ref G g)
	{
		scc = SCC(g);
		comp = Components(scc.condensate);
		top = TopOrder(scc.condensate);

		down.resize(scc.condensate.n);
		foreach(a; top.verts)
			foreach(b; scc.condensate.succ(a))
				down[b] = max(down[b], 1+down[a]);

		up.resize(scc.condensate.n);
		foreach_reverse(a; top.verts)
			foreach(b; scc.condensate.succ(a))
				up[a] = max(up[a], 1+up[b]);

		/*for(int a = 0; a < scc.condensate.n; ++a)
			foreach(b; scc.condensate.succ(a))
				assert(down[a] < down[b] && up[a] > up[b]);*/

		timeA.resize(scc.condensate.n);
		timeB.resize(scc.condensate.n);

		for(int i = 0; i < scc.condensate.n; ++i)
			if(down[i] == 0) // roots only
				dfs(i);
	}

	private void dfs(int a)
	{
		if(timeA[a] != 0)
			return;
		timeA[a] = time++;
		foreach(b; scc.condensate.succ(a))
			dfs(b);
		timeB[a] = time++;
	}

	int n() const pure nothrow @property
	{
		return cast(int)scc.comp.length;
	}

	/**
	 * tests if there is a path a->b.
	 * returns false if the answer can not be obtained in O(1).
	 */
	bool hasPathFast(int a, int b) const
	{
		return hasPathEx(a, b) == 1;
	}

	/**
	 * tests if there is a path a->b.
	 * 0 = no, 1 = yes, 2 = unknown
	 */
	int hasPathEx(int a, int b) const
	{
		// reduce to strongly connected components
		a = scc.comp[a];
		b = scc.comp[b];

		// inside the same SCC -> true
		if(a == b)
			return 1;

		// disconnected components -> false
		if(!comp.isConnected(a, b))
			return 0;

		// topological order wrong -> false
		if(!(top.order[a] < top.order[b]))
			return 0;

		// level order wrong -> false
		// NOTE: this is not always strictly stronger than top-order
		if(!(down[a] < down[b] && up[a] > up[b]))
			return 0;

		// part of the forest -> true
		if(timeA[a] < timeA[b] && timeA[b] < timeB[a])
			return 1;

		// unknown
		return 2;
	}

	/**
	 * number of true/false/unknown paths
	 * (debugging / performance measurements only)
	 */
	void testQuality(ref size_t no, ref size_t yes, ref size_t unknown)
	{
		no = yes = unknown = 0;

		for(int a = 0; a < n; ++a)
			for(int b = 0; b < n; ++b)
				final switch(hasPathEx(a, b))
				{
					case 0: no++; break;
					case 1: yes++; break;
					case 2: unknown++; break;
				}

		assert(yes + no + unknown == n*n);
	}
}

/**
 * Compute a topological order of the graph.
 * If it contains cycles, there is no topological order, but an approximation
 * of one is computed which still might be useful for heuristical algorithms.
 */
struct TopOrder
{
	Array!int order; // position of vertex in the order
	Array!int verts; // vertices in topological order

	@disable this();
	@disable this(this);

	this(G)(auto ref G g)
	{
		order.resize(g.n);
		verts.resize(g.n);

		auto v = BitArray(g.n);
		int cnt = g.n;

		void dfs(int a)
		{
			if(v[a])
				return;
			v[a] = true;
			foreach(e; g.succ(a))
				dfs(e);
			order[a] = --cnt;
			verts[order[a]] = a;
		}

		for(int a = 0; a < g.n; ++a)
			dfs(a);
		assert(cnt == 0);

		/*for(int a = 0; a < g.n; ++a)
			foreach(b; g.succ(a))
				assert(order[a] < order[b]);*/
	}
}

/** compute connected components of undirected graph */
struct Components
{
	Array!int comp;
	int nComps;

	this(G)(auto ref G g)
	{
		auto uf = UnionFind(g.n);
		for(int a = 0; a < g.n; ++a)
			foreach(b; g.succ(a))
				uf.join(a, b);
		comp = uf.components(nComps);
	}

	bool isConnected(int a, int b) const pure nothrow
	{
		return comp[a] == comp[b];
	}
}

/** compute strongly connected components of a directed graph */
struct SCC
{
	Array!int comp; // old vertex -> SCC
	Array!int compSize; // size of individual SCC's
	Graph!true condensate; // condensed graph, i.e. one vertex per SCC

	this(G)(auto ref G g)
	{
		auto visited = BitArray(g.n);
		auto back = Array!int(g.n);
		comp.resize(g.n);
		Array!int stack;
		int cnt = 0;

		void dfs(int v)
		{
			if(visited[v])
				return;
			visited[v] = true;

			int x = back[v] = cnt++;

			stack.pushBack(v);

			foreach(w; g.succ(v))
			{
				dfs(w);
				x = min(x, back[w]);
			}

			if(x < back[v])
			{
				back[v] = x;
				return;
			}

			int size = 0;
			while(true)
			{
				++size;
				int t = stack.popBack();
				back[t] = 999999999;
				comp[t] = nComps;

				if(t == v)
					break;
			}

			compSize.pushBack(size);
		}

		for(int i = 0; i < g.n; ++i)
			dfs(i);

		// build condensed graph
		condensate = Graph!true(nComps);
		for(int i = 0; i < g.n; ++i)
			foreach(e; g.succ(i))
				if(comp[i] != comp[e])
					condensate.addEdge(comp[i], comp[e]);
		condensate.reduce;
	}

	/** number of strongly connected components */
	int nComps() const nothrow @property
	{
		return cast(int)compSize.length;
	}

	/** size of the larges component */
	int largestComponentSize() const nothrow @property
	{
		int x = 0;
		foreach(s; compSize[])
			x = max(x, s);
		return x;
	}
}


//////////////////////////////////////////////////////////////////////
/// shortest path algorithms
//////////////////////////////////////////////////////////////////////

struct FloydWarshall(T)
{
	Array2!T dist;
	T diameter;

	@disable this();
	@disable this(this);

	this(G)(auto ref G g)
	{
		dist = Array2!T(g.n, g.n, T.max/2);
		for(int i = 0; i < g.n; ++i)
			foreach(j; g.succ(i))
				dist[i, j] = 1;

		for(int b = 0; b < g.n; ++b)
			for(int a = 0; a < g.n; ++a)
				if(dist[a,b] < T.max/2)
					for(int c = 0; c < g.n; ++c)
						dist[a,c] = min(dist[a,c], dist[a,b]+dist[b,c]);

		foreach(a, b, ref d; dist[])
			if(d == T.max/2)
				d = T.max;

		diameter = 0;
		foreach(a, b, d; dist[])
			diameter = max(diameter, d);
	}
}
