module algo.graph;

private import std.process : executeShell;
private import std.random;
private import std.stdio;
private import std.typetuple;
private import std.typecons;
private import std.math : sqrt, isNaN;
private import std.conv : to;
private import std.algorithm : min, max, move, sort;
private import jive.array;
private import jive.bitarray;

/**
 * directed or undirected (multi-)graph with optionally labeled edges
 * TODO: decide what type to use for vertex identifiers (int/long/size_t ?).
 */
class Graph(bool directed = false, Label...)
{
	static struct HalfEdge
	{
		int to;
		Label label;

		alias to this;
	}

	static struct Edge
	{
		int from;
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
	this(size_t n = 0)
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

	/** ignores edge labels */
	int diameter() const @property
	{
		return new FloydWarshall!int(this).diameter;
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
	bool hasEdge(int from, int to) const nothrow
	{
		foreach(e; succ(from))
			if(e.to == to)
				return true;
		return false;
	}

	/**
	 * tests if there is a path a->b.
	 * Uses a trivial O(n+m) dfs. For a more sophisticated
	 * approach, use the Reachability class
	 */
	bool hasPath(int from, int to) const
	{
		if(from == to)
			return true;

		auto v = BitArray(n);
		Array!int stack;
		stack.pushBack(from);
		v[from] = true;

		while(!stack.empty)
		{
			int a = stack.popBack;
			foreach(b; succ(a))
				if(b.to == to)
					return true;
				else if(!v[b.to])
				{
					v[b.to] = true;
					stack.pushBack(b.to);
				}
		}

		return false;
	}

	/** number of path-connected pairs of vertices. */
	long pathCount() const
	{
		size_t cnt = 0;
		BitArray v;

		for(int a = 0; a < n; ++a)
		{
			reachable(a, v);
			cnt += v.count(true) - 1; // "-1" for the a->a loop
		}

		return cnt;
	}

	/** determines the vertices reachable from a */
	BitArray reachable(int a) const
	{
		BitArray v;
		reachable(a, v);
		return v;
	}

	/** ditto (can be used to prevent heap alocation) */
	void reachable(int a, ref BitArray v) const
	{
		v.resize(n);
		v.reset();
		Array!int stack;

		stack.pushBack(a);
		v[a] = true;

		while(!stack.empty)
		{
			int b = stack.popBack;
			foreach(c; succ(b))
				if(!v[c])
				{
					v[c] = true;
					stack.pushBack(c);
				}
		}
	}

	/** range of successors of a vertex */
	const(HalfEdge)[] succ(int x) const nothrow
	{
		return g[x][];
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
			foreach(edge; succ(a))
				if(directed || a <= edge.to)
					static if(Label.length)
						f.writefln("%s%s%s [label=\"%s\"]", a, directed?"->":"--", edge.to, to!string(edge.label));
					else
						f.writefln("%s%s%s", a, directed?"->":"--", edge.to);
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

Graph!dir completeGraph(bool dir = false)(int n)
{
	auto g = new Graph!dir(n);
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
	auto g = new Graph!dir(n);
	for(int i = 1; i < n; ++i)
		g.addEdge(i-1, i);
	return g;
}

Graph!dir randomGraph(bool dir = false)(int n, float ratio)
{
	if((dir && ratio >= n-1) || (!dir && ratio >= (n-1)/2))
		return completeGraph!dir(n);

	auto g = new Graph!dir(n);
	while(g.ratio < ratio)
	{
		int a = uniform(0, n);
		int b = uniform(0, n);
		if(a != b && !g.hasEdge(a, b))
			g.addEdge(a, b);
	}
	return g;
}

Graph!dir randomTree(bool dir = false)(int n)
{
	auto g = new Graph!dir(n);
	for(int i = 1; i < n; ++i)
		g.addEdge(i, uniform(0, i));
	return g;
}

Graph!dir binaryLeftTree(bool dir = false)(int n)
{
	auto g = new Graph!dir(n);
	for(int i = 1; i < n; ++i)
		g.addEdge(i, (i-1)/2);
	return g;
}

Graph!(false, float) randomEuclideanGraph(int n, float size, float cutoff)
{
	auto g = new Graph!(false, float)(n);
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

class FloydWarshall(T)
{
	Array2!T dist;
	T diameter;

	this(G)(G g)
	{
		dist = Array2!T(g.n, g.n, T.max/2);
		for(int i = 0; i < g.n; ++i)
			foreach(e; g.succ(i))
				dist[i, e.to] = 1;

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

/**
 * Compute a topological order of the graph.
 * If it contains cycles, there is no topological order, but an approximation
 * of one is computed which still might be useful for heuristical algorithms.
 */
class TopOrder
{
	Array!int order; // position of vertex in the order
	Array!int verts; // vertices in topological order

	this(G)(G g)
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
				dfs(e.to);
			order[a] = --cnt;
			verts[order[a]] = a;
		}

		for(int a = 0; a < g.n; ++a)
			dfs(a);
		assert(cnt == 0);
	}
}

/** compute strongly connected components of a directed graph */
class SCC : Graph!true
{
	Array!int comp; // old vertex -> SCC
	Array!int compSize; // size of individual SCC's

	int nComps() const nothrow @property
	{
		return cast(int)compSize.length;
	}

	this(G)(G g)
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
				dfs(w.to);
				x = min(x, back[w.to]);
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
		super(nComps);
		for(int i = 0; i < g.n; ++i)
			foreach(e; g.succ(i))
				if(comp[i] != comp[e.to])
					addEdge(comp[i], comp[e.to]);
		this.reduce;
	}
}

/**
 * Answers reachability queries for directed graphs.
 *
 * Inteded for large graphs where building the complete transitive closure is
 * not an option. Non-exact queries may return false even if a path exists, but
 * only ever take O(1) time. Positive answers are always correct.
 */
class Reachability
{
	SCC s;
	Array!long timeA;
	Array!long timeB;
	Array!int nodes;
	TopOrder top;

	this(G)(G g)
	{
		s = new SCC(g);
		top = new TopOrder(s);

		timeA.resize(s.n, -1);
		timeB.resize(s.n, -1);

		void dfs(int a)
		{
			if(timeA[a] != -1)
				return;
			timeA[a] = nodes.length;
			nodes.pushBack(a);

			foreach(b; s.succ(a))
				dfs(b);
			timeB[a] = nodes.length;
		}

		// sort edges topological (improvement of this is quite small on random graphs)
		for(int i = 0; i < s.n; ++i)
			sort!((a,b)=> top.order[a.to] < top.order[b.to])(s.g[i][]);

		// build tree in topological order
		foreach(i; top.verts)
			dfs(i);
	}

	/** checks wether there is a path a -> b */
	bool hasPath(bool exact)(int a, int b) const
	{
		// replace vertices by there SCC's
		a = s.comp[a];
		b = s.comp[b];

		// part of thre forest -> true (including the case a == b)
		if(timeA[a] <= timeA[b] && timeA[b] <= timeB[a])
			return true;

		// contradicted by topological order -> false
		if(top.order(a) > top.order(b))
			return false;

		// answer not clear -> revert to actual dfs (on the SCC graph)
		static if(exact)
			return s.hasPath(a, b);
		else
			return false;
	}

	/** count how many pairs (a,b) are connected by a path */
	long pathCount(bool exact)() const @property
		if(exact == true)
	{
		long cnt = 0;
		BitArray v;
		for(int i = 0; i < s.n; ++i)
		{
			s.reachable(i, v);
			for(int j = 0; j < s.n; ++j)
				if(v[j])
					cnt += s.compSize[i] * s.compSize[j];
		}
		cnt -= s.comp.length;
		return cnt;
	}

	long pathCount(bool exact)() const @property
		if(exact == false)
	{
		long cnt = 0;
		for(int i = 0; i < s.n; ++i)
			foreach(j; nodes[timeA[i]..timeB[i]])
				cnt += s.compSize[i] * s.compSize[j];
		cnt -= s.comp.length;
		return cnt;
	}
}
