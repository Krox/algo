module algo.closestpair;

private import std.algorithm;
private import std.math : sqrt;
private import jive.set;

import jive.array : Array;
import algo.geometry;

/**
 * Returns the closest pair of points. The problem is not that interesting,
 * but the algorithm is kinda nice as it achieves O(n log n) runtime.
 */
Tuple!(Vec2!T,Vec2!T) closestPair(T)(Array!(Vec2!T) points)
	if(is(T == double)) // TODO: figure out how to do it for exact types
{
	if(points.length < 2)
		throw new Exception("need at least two points for a pair");

	sort!"a.x < b.x"(points[]);

	Set!(Vec2!T, "a.y < b.y || (a.y == b.y && a.x < b.x)") sss;

	size_t tail = 0; // next one to be deleted

	T best = T.max/2; // rounded up
	Tuple!(Vec2!T,Vec2!T) r;

	foreach(p; points[])
	{
		while(p.x - points[tail].x > best)
			sss.remove(points[tail++]);
		foreach(q; sss.range!"[]"(p-Vec2!T(best,best),p+Vec2!T(best,best)))
			if(sqrt(dot!(T,2)(p-q, p-q)) < best)
			{
				best = sqrt(dot!(T,2)(p-q, p-q));
				r[0] = q;
				r[1] = p;
			}
		sss.add(p);
	}

	return r;
}
