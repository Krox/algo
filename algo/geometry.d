module algo.geometry;


private import std.traits;
private import std.conv;

/**
 * Matrix/Vector of (small) fixed size.
 * Mainly intended for geometry in 2D,3D,4D.
 */
struct Mat(T, size_t N, size_t M)
{
	T[N][M] c;

	enum Mat zero = constant(0);
	enum Mat identity = _identity();

	static Mat constant(T v) pure
	{
		Mat r;
		foreach(n; 0..N)
			foreach(m; 0..M)
				r[n,m] = v;
		return r;
	}

	static Mat _identity() pure
	{
		Mat r;
		foreach(n; 0..N)
			foreach(m; 0..M)
				r[n,m] = n==m ? 1 : 0;
		return r;
	}

	ref T opIndex(size_t n, size_t m)
	{
		return c[m][n];
	}

	ref const(T) opIndex(size_t n, size_t m) const
	{
		return c[m][n];
	}

	static if(M == 1)
	{
		ref inout(T) opIndex(size_t n) inout
		{
			return c[0][n];
		}
	}

	static if(N == 1)
	{
		ref inout(T) opIndex(size_t m) inout
		{
			return c[m][0];
		}
	}

	static if(N == 1 || M == 1)
	{
		ref inout(T) x() inout @property { return this[0]; }
		ref inout(T) y() inout @property { return this[1]; }
		ref inout(T) z() inout @property { return this[2]; }
		ref inout(T) w() inout @property { return this[3]; }

		this(T x) { this.x = x; }
		this(T x, T y) { this.x = x; this.y = y; }
		this(T x, T y, T z) { this.x = x; this.y = y; this.z = z; }
		this(T x, T y, T z, T w) { this.x = x; this.y = y; this.z = z; this.w = w; }
	}

	void opAddAssign(const Mat b)
	{
		foreach(n; 0..N)
			foreach(m; 0..M)
				 this[n,m] += b[n,m];
	}

	void opSubAssign(const Mat b)
	{
		foreach(n; 0..N)
			foreach(m; 0..M)
				 this[n,m] -= b[n,m];
	}

	void opMulAssign(T s)
	{
		foreach(n; 0..N)
			foreach(m; 0..M)
				 this[n,m] *= s;
	}

	Mat opMul(T s) const
	{
		Mat r = this;
		r *= s;
		return r;
	}

	static if(isFloatingPoint!T)
	{
		void opDivAssign(T s)
		{
			opMulAssign(1.0/s);
		}

		Mat opDiv(T s) const
		{
			return opMul(1.0/s);
		}
	}

	Mat opAdd(const Mat b) const
	{
		Mat r = this;
		r += b;
		return r;
	}

	Mat opSub(const Mat b) const
	{
		Mat r = this;
		r -= b;
		return r;
	}

	Mat!(T,N,K) opMul(size_t K)(ref const Mat!(T,M,K) rhs) const
	{
		auto r = Mat!(T,N,K).zero;
		foreach(n; 0..N)
			foreach(k; 0..K)
				foreach(m; 0..M)
					r[n,k] += this[n,m]*rhs[m,k];
		return r;
	}

	static if(N==M) Mat pow(long exp) const
	{
		assert(exp >= 0);
		Mat r = identity;
		Mat base = this;

		while(exp)
		{
			if(exp & 1)
				r = r * base;
			exp >>= 1;
			base = base * base;
		}
		return r;
	}

	string toString() const @property
	{
		return to!string(c);	// TODO: transpose this
	}
}

alias Vec(T, size_t N)  = Mat!(T, N, 1);

alias Mat2(T) = Mat!(T, 2, 2);
alias Mat3(T) = Mat!(T, 2, 2);
alias Mat4(T) = Mat!(T, 2, 2);

alias Vec2(T) = Vec!(T, 2);
alias Vec3(T) = Vec!(T, 3);
alias Vec4(T) = Vec!(T, 4);

/**
 * (standard/euclidean) inner product of two vectors
 */
T dot(T, size_t N)(Vec!(T,N) a, Vec!(T,N) b)
{
	T sum = T(0);
	for(size_t i = 0; i < N; ++i)
		sum = sum + a[i]*b[i];
	return sum;
}

/**
 * (z-component of the) crossproduct of two vectors
 */
T cross(T)(Vec2!T a, Vec2!T b)
{
	return a[0] * b[1] - a[1] * b[0];
}

/**
 * return true if (a,b,c) is a non-degenerate triangle in counter-clockwise orientation
 */
bool ccw(T)(Vec2!T a, Vec2!T b, Vec2!T c)
{
	return cross!T(b-a, c-a) > 0;
}
