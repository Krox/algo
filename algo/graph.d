module algo.graph;

private import std.process : system;
private import std.random;
private import std.stdio;
private import std.typetuple;
private import std.typecons;
private import std.math : sqrt, isnan;
private import std.conv : to;

import jive.array;

/**
 * directed or undirected (multi-)graph with optionally labeled edges
 */
class Graph(bool directed, Label...)
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
	size_t numVertices() const @property
	{
		return g.length;
	}

	/** number of edges */
	size_t numEdges() const @property
	{
		return edgeCount;
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

	/** range of successors of a vertex */
	const(HalfEdge)[] succ(int x) const
	{
		return g[x][];
	}

	/** write graph into a .dot file */
	void writeDot(string filename)
	{
		auto f = File(filename, "w");
		f.writefln("%s %s {", directed?"digraph":"graph", name);

		//f.writefln("node [label=\"\\N\"];");
		if(!isnan(sizeX) && !isnan(sizeY))
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
		system(layout~" "~"tmp_"~name~".dot -Tsvg > "~"tmp_"~name~".svg");
		system("eog "~"tmp_"~name~".svg");
		system("rm "~"tmp_"~name~".dot "~"tmp_"~name~".svg");
	}

	static if(Label.length == 0)
	{
		static Graph createLine(int n)
		{
			auto g = new Graph;
			auto a = g.addVertex();
			for(int i = 1; i < n; ++i)
			{
				auto b = g.addVertex();
				g.addEdge(a,b);
				a = b;
			}
			return g;
		}

		static Graph createRandomTree(int n)
		{
			auto g = new Graph;
			g.addVertex();
			for(int i = 1; i < n; ++i)
			{
				auto a = g.addVertex();
				g.addEdge(uniform(0,i), a);
			}
			return g;
		}

		static Graph createLeftBinaryTree(int n)
		{
			auto g = new Graph;
			g.addVertex();
			for(int i = 1; i < n; ++i)
			{
				auto a = g.addVertex();
				g.addEdge((i-1)/2, a);
			}
			return g;
		}
	}

	static if(Label.length == 1 && is(Label[0] == float))
	{
		static Graph createRandomEuclidean(int n, float size, float cutoff)
		{
			auto g = new Graph(n);
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
	}
}
