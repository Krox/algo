module algo.convex;

private import std.algorithm;

import jive.array;
import algo.geometry;


/**
 * Convex hull of a set of points, implemented using the monotone chains
 * algorithm. Output is sorted ccw, starting with the left(/bottom) most
 * point. Returns the minimal convex hull, i.e. duplicate points and
 * colinear points on the hull are ignored.
 */
Array!(Vec2!T) convexHull(T)(Array!(Vec2!T) points)
{
	// sort points and remove duplicates
	sort!"a.x < b.x || (a.x == b.x && a.y < b.y)"(points[]);

	Array!(Vec2!T) r;

	// stupid special cases
	if(points.empty)
		return r;
	if(points[0] == points[$-1])
	{
		r.pushBack(points[0]);
		return r;
	}

	// lower half of hull
	foreach(p; points[])
	{
		while(r.length >= 2 && !ccw!T(r[$-2], r[$-1], p))
			r.popBack();
		r.pushBack(p);
	}
	r.popBack();

	auto low = r.length;

	// upper half of hull
	foreach_reverse(p; points[])
	{
		while(r.length-low >= 2 && !ccw!T(r[$-2], r[$-1], p))
			r.popBack();
		r.pushBack(p);
	}
	r.popBack();

	return r;
}
