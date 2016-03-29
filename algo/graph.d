module algo.graph;

private import std.process : executeShell;
private import std.random;
private import std.stdio;
private import std.typetuple;
private import std.typecons;
private import std.math : sqrt, isNaN;
private import std.conv : to;
private import std.algorithm : min, max;
private import jive.array;

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
	this(size_t numVertices = 0)
	{
		assert(numVertices <= int.max);
		g.resize(numVertices);
	}

	/** number of vertices */
	size_t numVertices() const nothrow @safe @property
	{
		return g.length;
	}

	/** number of edges */
	size_t numEdges() const nothrow @safe @property
	{
		return edgeCount;
	}

	/** maximum number of edges possible with current vertex count */
	size_t maxEdges() const nothrow @safe @property
	{
		static if(directed)
			return numVertices * (numVertices - 1);
		else
			return numVertices * (numVertices - 1) / 2;
	}

	/** numEdges / maxEdges */
	float density() const nothrow @safe @property
	{
		return cast(float)numEdges / maxEdges;
	}

	/** numEdges / numVertices */
	float ratio() const nothrow @safe @property
	{
		return cast(float)numEdges / numVertices;
	}

	/** ignores edge labels */
	int diameter() const @property
	{
		return new FloydWarshall!int(this).diameter;
	}

	/** add a new vertex, returns its number */
	int addVertex()
	{
		assert(numEdges < int.max);
		g.pushBack((Array!HalfEdge).init);
		return cast(int)(numVertices-1);
	}

	/** add a new edge */
	void addEdge(int from, int to, Label label)
	{
		assert(0 <= from && from < numVertices);
		assert(0 <= to && to < numVertices);
		g[from].pushBack(HalfEdge(to, label));
		if(!directed && to != from)
			g[to].pushBack(HalfEdge(from, label));
		++edgeCount;
	}

	/** test wether a certainedge exists */
	bool hasEdge(int from, int to) const nothrow
	{
		foreach(e; succ(from))
			if(e.to == to)
				return true;
		return false;
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

		for(int a = 0; a < numVertices; ++a)
			if(a < pos.length)
				f.writefln("%s [pos=\"%s,%s\"];", a, pos[a][0], pos[a][1]);
			else
				f.writefln("%s", a);

		for(int a = 0; a < numVertices; ++a)
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
		return completeGraph(n);

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
		dist = Array2!T(g.numVertices, g.numVertices, T.max/2);
		for(int i = 0; i < g.numVertices; ++i)
			foreach(e; g.succ(i))
				dist[i, e.to] = 1;

		for(int b = 0; b < g.numVertices; ++b)
			for(int a = 0; a < g.numVertices; ++a)
				if(dist[a,b] < T.max/2)
					for(int c = 0; c < g.numVertices; ++c)
						dist[a,c] = min(dist[a,c], dist[a,b]+dist[b,c]);

		foreach(a, b, ref d; dist[])
			if(d == T.max/2)
				d = T.max;

		diameter = 0;
		foreach(a, b, d; dist[])
			diameter = max(diameter, d);
	}
}
