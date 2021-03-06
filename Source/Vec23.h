/* 
 * heliangliang, USTB, 2012.03.05, CopyRight Reserved
*/

#ifndef VEC23_H_
#define VEC23_H_

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <cmath>
#include <iostream>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

template<typename _T>
class mat33;
typedef mat33<double>	mat33d;
typedef mat33<float>	mat33f;


template<typename _T>
class vec3;

typedef vec3<float>		vec3f;
typedef vec3<double>	vec3d;
typedef vec3<int>		vec3i;

template<typename _T>
class vec2;

typedef vec2<float>		vec2f;
typedef vec2<double>	vec2d;
typedef vec2<int>		vec2i;


template<int _D, typename _T>
class vec_dim {
public:
	typedef vec3<_T> T;
};

template<typename _T>
class vec_dim<2, _T> {
public:
	typedef vec2<_T> T;
};

//3x3 matrix
template<int _D,typename _T>
class mat_dim{
public: 
	typedef mat33<_T> T;
};


template<typename _T>
class vec3 {
public:
	typedef _T type;
	_T x, y, z;
	// Constructor,Destructor
	//inli ne vec3(): x(_T(0)), y(_T(0)), z(_T(0)) { }
    inline vec3() { }
	inline explicit vec3( _T a ): x(a), y(a), z(a) { }
	inline vec3( _T a, _T b, _T c ): x(a), y(b), z(c) { }
	inline vec3( const vec3<_T>& a ): x(a.x), y(a.y), z(a.z) { }
	inline explicit vec3( const vec2<_T>& a ): x(a.x), y(a.y), z(0) { }
	inline ~vec3() { }

	// ==, <, >, <=, >=
	inline bool operator==(const vec3<_T>& a) const
		{ return (x==a.x && y==a.y && z==a.z); }
    inline bool operator<(const vec3<_T>& a) const
        { return (x<a.x && y<a.y && z<a.z); }
    inline bool operator>(const vec3<_T>& a) const
        { return (x>a.x && y>a.y && z>a.z); }
    inline bool operator<=(const vec3<_T>& a) const
        { return (x<=a.x && y<=a.y && z<=a.z); }
    inline bool operator>=(const vec3<_T>& a) const
        { return (x>=a.x && y>=a.y && z>=a.z);; }
	// Set
	inline void set( _T a )
		{ x=a; y=a; z=a; }
	inline void set( _T a, _T b, _T c )
		{ x=a; y=b; z=c; }
	// []
	inline _T& operator[]( int i )
		{ return *(&x+i); }
	inline const _T& operator[]( int i ) const
		{ return *(&x+i); }
	// -
	inline vec3<_T> operator-() const
		{ return vec3<_T>(-x,-y,-z); }
	// =
	inline const vec3<_T>& operator=(const vec3<_T>& a)
		{ /*if(this==&a)return *this;*/ x=a.x; y=a.y; z=a.z; return *this; }
	// +=, +
	inline const vec3<_T>& operator+=( _T a )
		{ x+=a; y+=a; z+=a; return *this; }
	inline const vec3<_T>& operator+=( const vec3<_T>& a )
		{ x+=a.x; y+=a.y; z+=a.z; return *this; }
	inline vec3<_T> operator+( _T a ) const
		{ return vec3<_T>(x+a,y+a,z+a); }
	inline vec3<_T> operator+( const vec3<_T>&a ) const
		{ return vec3<_T>(x+a.x,y+a.y,z+a.z); }
	// -=, -
	inline const vec3<_T>& operator-=( _T a )
		{ x-=a; y-=a; z-=a; return *this; }
	inline const vec3<_T>& operator-=( const vec3<_T>& a )
		{ x-=a.x; y-=a.y; z-=a.z; return *this; }
	inline vec3<_T> operator-( _T a ) const
		{ return vec3<_T>(x-a,y-a,z-a); }
	inline vec3<_T> operator-( const vec3<_T>&a ) const
		{ return vec3<_T>(x-a.x,y-a.y,z-a.z); }
	// *=, *
	inline const vec3<_T>& operator*=( _T a )
		{ x*=a; y*=a; z*=a; return *this; }
	inline const vec3<_T>& operator*=( const vec3<_T>& a )
		{ x*=a.x; y*=a.y; z*=a.z; return *this; }
	inline vec3<_T> operator*( _T a ) const
		{ return vec3<_T>(x*a,y*a,z*a); }
	inline vec3<_T> operator*( const vec3<_T>&a ) const
		{ return vec3<_T>(x*a.x,y*a.y,z*a.z); }
	// /=, /
	inline const vec3<_T>& operator/=( _T a )
		{ x/=a; y/=a; z/=a; return *this; }
	inline const vec3<_T>& operator/=( const vec3<_T>& a )
		{ x/=a.x; y/=a.y; z/=a.z; return *this; }
	inline vec3<_T> operator/( _T a ) const
		{ return vec3<_T>(x/a,y/a,z/a); }
	inline vec3<_T> operator/( const vec3<_T>&a ) const
		{ return vec3<_T>(x/a.x,y/a.y,z/a.z); }

	// Dot product
	inline _T dot( const vec3<_T>&a ) const
		{ return (x*a.x+y*a.y+z*a.z); }
	// Cross product
	inline vec3<_T> cross( const vec3<_T>&a ) const
		{ return vec3<_T>( y*a.z-z*a.y, z*a.x-x*a.z, x*a.y-y*a.x ); }
	inline vec3<_T> cross( const vec2<_T>&a ) const
		{ return vec3<_T>( -z*a.y, z*a.x, x*a.y-y*a.x ); }




	// Product of a's transpose and t
	inline mat33<_T> pro_transpose(const vec3<_T>& a) const;
	inline mat33<_T> pro_trans_self() const;
	// vec3 t; mat33 a; a * t
	inline vec3<_T> pro_matrix_left(const mat33<_T>& a) const;
	// vec3 t; mat33 a; t^T * a
	inline vec3<_T> pro_matrix_right(const mat33<_T>& a) const;
	


	// Normalize
	inline const vec3<_T>& normalize()
		{ _T l=x*x+y*y+z*z; if(l>_T(0)){l=std::sqrt(l);*this/=l;} return *this; }
	// Length,矩阵的模
	inline _T length() const
		{ return std::sqrt(len_square()); }
	inline _T len_square() const
		{ return x*x+y*y+z*z; }
	inline _T distance(const vec3<_T>& a) const
		{ return std::sqrt(dis_square(a)); }
	inline _T dis_square(const vec3<_T>& a) const
		{ return (x-a.x)*(x-a.x)+(y-a.y)*(y-a.y)+(z-a.z)*(z-a.z); }

	// Pointer
	inline _T* ptr()
		{ return &x; }
	inline const _T* ptr() const
		{ return &x; }

	// Max, Min
	inline _T max() const
		{ return std::max(x,std::max(y,z)); }
	inline _T min() const
		{ return std::min(x,std::min(y,z)); }
	inline vec3<_T> max(const vec3<_T>&a) const
		{ return vec3<_T>(std::max(x,a.x),std::max(y,a.y),std::max(z,a.z)); }
	inline vec3<_T> min(const vec3<_T>&a) const
		{ return vec3<_T>(std::min(x,a.x),std::min(y,a.y),std::min(z,a.z)); }

    // Volume
    inline _T volume() const
        { return x*y*z; }
	// inside
	inline bool inside(const vec3<_T>& pmin, const vec3<_T>& pmax) const
		{ return pmin<=*this && *this<=pmax; }

	// To other type
    template<typename _toT>
    inline vec3<_toT> to() const
        { return vec3<_toT>( _toT(x), _toT(y), _toT(z) ); }
    /*inline vec3<int> to_int() const
    { return vec3<int>( int(x), int(y), int(z) ); }
    inline vec3<float> to_float() const
    { return vec3<float>( float(x), float(y), float(z) ); }
    inline vec3<double> to_double() const
    { return vec3<double>( double(x), double(y), double(z) ); }*/

	// << and >>
	friend inline std::ostream& operator<<( std::ostream& os, const vec3<_T>& a )
		{ os << a.x << ' ' << a.y << ' ' << a.z; return os; }
	friend inline std::istream& operator>>( std::istream& is, vec3<_T>& a )
		{ is >> a.x >> a.y >> a.z; return is; }
	friend inline std::ofstream& operator<<( std::ofstream& os, const vec3<_T>& a )
		{ os << a.x << ' ' << a.y << ' ' << a.z; return os; }
	friend inline std::ifstream& operator>>( std::ifstream& is, vec3<_T>& a )
		{ is >> a.x >> a.y >> a.z; return is; }

	// Data
	static const vec3<_T> O/*, E, Ex, Ey, Ez*/;
	static const int dim = 3;

};

//template<typename _T>
//const vec3<_T> vec3<_T>::E=vec3<_T>(_T(1));
//template<typename _T>
//const vec3<_T> vec3<_T>::Ex=vec3<_T>(_T(1),_T(0),_T(0));
//template<typename _T>
//const vec3<_T> vec3<_T>::Ey=vec3<_T>(_T(0),_T(1),_T(0));
//template<typename _T>
//const vec3<_T> vec3<_T>::Ez=vec3<_T>(_T(0),_T(0),_T(1));
template<typename _T>
const vec3<_T> vec3<_T>::O=vec3<_T>(_T(0));



template<typename _T>
class vec2 {
public:
	typedef _T type;
	_T x, y;
	// Constructor,Destructor
	inline vec2() { }
	inline explicit vec2( _T a ): x(a), y(a) { }
	inline vec2( _T a, _T b ): x(a), y(b) { }
	inline vec2( const vec2<_T>& a ): x(a.x), y(a.y) { }
	inline explicit vec2( const vec3<_T>& a ): x(a.x), y(a.y) { }
	inline ~vec2() { }
	// ==, <, >, <=, >=
	inline bool operator==(const vec2<_T>& a) const
		{ return (x==a.x && y==a.y); }
	inline bool operator<(const vec2<_T>& a) const
		{ return (x<a.x && y<a.y); }
	inline bool operator>(const vec2<_T>& a) const
		{ return (x>a.x && y>a.y); }
	inline bool operator<=(const vec2<_T>& a) const
		{ return (x<=a.x && y<=a.y); }
	inline bool operator>=(const vec2<_T>& a) const
		{ return (x>=a.x && y>=a.y); }
	// Set
	inline void set( _T a )
		{ x=a; y=a; }
	inline void set( _T a, _T b )
		{ x=a; y=b; }
	// []
	inline _T& operator[]( int i )
		{ return *(&x+i); }
	inline const _T& operator[]( int i ) const
		{ return *(&x+i); }
	// -
	inline vec2<_T> operator-() const
		{ return vec2<_T>(-x,-y); }
	// =
	inline const vec2<_T>& operator=(const vec2<_T>& a)
		{ /*if(this==&a)return *this;*/ x=a.x; y=a.y; return *this; }
	// +=, +
	inline const vec2<_T>& operator+=( _T a )
		{ x+=a; y+=a; return *this; }
	inline const vec2<_T>& operator+=( const vec2<_T>& a )
		{ x+=a.x; y+=a.y; return *this; }
	inline vec2<_T> operator+( _T a ) const
		{ return vec2<_T>(x+a,y+a); }
	inline vec2<_T> operator+( const vec2<_T>&a ) const
		{ return vec2<_T>(x+a.x,y+a.y); }
	// -=, -
	inline const vec2<_T>& operator-=( _T a )
		{ x-=a; y-=a; return *this; }
	inline const vec2<_T>& operator-=( const vec2<_T>& a )
		{ x-=a.x; y-=a.y; return *this; }
	inline vec2<_T> operator-( _T a ) const
		{ return vec2<_T>(x-a,y-a); }
	inline vec2<_T> operator-( const vec2<_T>&a ) const
		{ return vec2<_T>(x-a.x,y-a.y); }
	// *=, *
	inline const vec2<_T>& operator*=( _T a )
		{ x*=a; y*=a; return *this; }
	inline const vec2<_T>& operator*=( const vec2<_T>& a )
		{ x*=a.x; y*=a.y; return *this; }
	inline vec2<_T> operator*( _T a ) const
		{ return vec2<_T>(x*a,y*a); }
	inline vec2<_T> operator*( const vec2<_T>&a ) const
		{ return vec2<_T>(x*a.x,y*a.y); }
	// /=, /
	inline const vec2<_T>& operator/=( _T a )
		{ x/=a; y/=a; return *this; }
	inline const vec2<_T>& operator/=( const vec2<_T>& a )
		{ x/=a.x; y/=a.y; return *this; }
	inline vec2<_T> operator/( _T a ) const
		{ return vec2<_T>(x/a,y/a); }
	inline vec2<_T> operator/( const vec2<_T>&a ) const
		{ return vec2<_T>(x/a.x,y/a.y); }

	// Dot product
	inline _T dot( const vec2<_T>&a ) const
		{ return (x*a.x+y*a.y); }
	// Cross product, return vec3 type
	//inline vec2<_T> cross( const vec2<_T>&a ) const
	//	{ return vec2<_T>( 0 ); }
	inline vec3<_T> cross( const vec2<_T>&a ) const
		{ return vec3<_T>( 0, 0, x*a.y-y*a.x ); }

	// Normalize
	inline const vec2<_T>& normalize()
		{ _T l=x*x+y*y; if(l>_T(0)){l=std::sqrt(l);*this/=l;} return *this; }
	// Length
	inline _T length() const
		{ return std::sqrt(len_square()); }
	inline _T len_square() const
		{ return x*x+y*y; }
	inline _T distance(const vec2<_T>& a) const
		{ return std::sqrt(dis_square(a)); }
	inline _T dis_square(const vec2<_T>& a) const
		{ return (x-a.x)*(x-a.x)+(y-a.y)*(y-a.y); }

	// Pointer
	inline _T* ptr()
		{ return &x; }
	inline const _T* ptr() const
		{ return &x; }
	// Max, Min
	inline _T max() const
		{ return std::max(x,y); }
	inline _T min() const
		{ return std::min(x,y); }
	inline vec2<_T> max(const vec2<_T>&a) const
		{ return vec2<_T>(std::max(x,a.x),std::max(y,a.y)); }
	inline vec2<_T> min(const vec2<_T>&a) const
		{ return vec2<_T>(std::min(x,a.x),std::min(y,a.y)); }
	// Volume
	inline _T volume() const
		{ return x*y; }
	// inside
	inline bool inside(const vec2<_T>& pmin, const vec2<_T>& pmax) const
		{ return pmin<=*this && *this<=pmax; }

	// To other type
	template<typename _toT>
	inline vec2<_toT> to() const
		{ return vec2<_toT>( _toT(x), _toT(y) ); }

	// << and >>
	friend inline std::ostream& operator<<( std::ostream& os, const vec2<_T>& a )
		{ os << a.x << ' ' << a.y; return os; }
	friend inline std::istream& operator>>( std::istream& is, vec2<_T>& a )
		{ is >> a.x >> a.y; return is; }
	friend inline std::ofstream& operator<<( std::ofstream& os, const vec2<_T>& a )
		{ os << a.x << ' ' << a.y; return os; }
	friend inline std::ifstream& operator>>( std::ifstream& is, vec2<_T>& a )
		{ is >> a.x >> a.y; return is; }

	// Data
	static const vec2<_T> O/*, E, Ex, Ey*/;
	static const int dim = 2;

};

//template<typename _T>
//const vec2<_T> vec2<_T>::E=vec2<_T>(_T(1));
//template<typename _T>
//const vec2<_T> vec2<_T>::Ex=vec2<_T>(_T(1),_T(0));
//template<typename _T>
//const vec2<_T> vec2<_T>::Ey=vec2<_T>(_T(0),_T(1));
template<typename _T>
const vec2<_T> vec2<_T>::O=vec2<_T>(_T(0));



template<typename _T>
class mat33 {
public:
	typedef _T type;
	vec3<_T> U, V, W;
	// Constructor,Destructor
	//inline mat33(): U(_T(0)), V(_T(0)), W(_T(0)) { }
    inline mat33() { } 
	explicit inline mat33( _T a ): U(a), V(a), W(a) { }
	inline mat33( _T a0, _T a1, _T a2, _T a3, _T a4, _T a5, _T a6, _T a7, _T a8 ):
		U(a0,a1,a2),V(a3,a4,a5),W(a6,a7,a8) { }
	explicit inline mat33( const vec3<_T>& a ): U(a), V(a), W(a) { }
	inline mat33( const vec3<_T>& a, vec3<_T>& b, vec3<_T>& c ): U(a), V(b), W(c) { }
	inline mat33( const mat33<_T>& a ): U(a.U), V(a.V), W(a.W) { }
	inline ~mat33() { }
	// Set
	inline void set( _T a )
		{ U.set(a); V.set(a); W.set(a); }
	inline void set( vec3<_T>& a )
		{ U=a; V=a; W=a; }
	inline void set( vec3<_T>& a, vec3<_T>& b, vec3<_T>& c )
		{ U=a; V=b; W=c; }
	// []
	inline vec3<_T>& operator[]( int i )
		{ return *(&U+i); }
	inline const vec3<_T>& operator[]( int i ) const
		{ return *(&U+i); }
	// -
	inline mat33<_T> operator-() const
		{ return mat33<_T>(-U,-V,-W); }
	// =
	inline const mat33<_T>& operator=(const mat33<_T>& a)
		{ if(this==&a)return *this; U=a.U; V=a.V; W=a.W; return *this; }
	// +=, +
	inline const mat33<_T>& operator+=( _T a )
		{ U+=a; V+=a; W+=a; return *this; }
	inline const mat33<_T>& operator+=( const mat33<_T>& a )
		{ U+=a.U; V+=a.V; W+=a.W; return *this; }
	inline mat33<_T> operator+( _T a ) const
		{ return mat33<_T>(U+a,V+a,W+a); }
	inline mat33<_T> operator+( const mat33<_T>&a ) const
		{ return mat33<_T>(U+a.U,V+a.V,W+a.W); }
	// -=, -
	inline const mat33<_T>& operator-=( _T a )
		{ U-=a; V-=a; W-=a; return *this; }
	inline const mat33<_T>& operator-=( const mat33<_T>& a )
		{ U-=a.U; V-=a.V; W-=a.W; return *this; }
	inline mat33<_T> operator-( _T a ) const
		{ return mat33<_T>(U-a,V-a,W-a); }
	inline mat33<_T> operator-( const mat33<_T>&a ) const
		{ return mat33<_T>(U-a.U,V-a.V,W-a.W); }
	// *=, *
	inline const mat33<_T>& operator*=( _T a )
		{ U*=a; V*=a; W*=a; return *this; }
	inline const mat33<_T>& operator*=( const mat33<_T>& a )
		{ U*=a.U; V*=a.V; W*=a.W; return *this; }
	inline mat33<_T> operator*( _T a ) const
		{ return mat33<_T>(U*a,V*a,W*a); }
	inline mat33<_T> operator*( const mat33<_T>&a ) const
		{ return mat33<_T>(U*a.U,V*a.V,W*a.W); }
	// /=, /
	inline const mat33<_T>& operator/=( _T a )
		{ U/=a; V/=a; W/=a; return *this; }
	inline const mat33<_T>& operator/=( const mat33<_T>& a )
		{ U/=a.U; V/=a.V; W/=a.W; return *this; }
	inline mat33<_T> operator/( _T a ) const
		{ return mat33<_T>(U/a,V/a,W/a); }
	inline mat33<_T> operator/( const mat33<_T>&a ) const
		{ return mat33<_T>(U/a.U,V/a.V,W/a.W); }
	
	// Product
	inline mat33<_T> product( const mat33<_T>& a ) const
		{ return mat33<_T>(
			U.x*a.U.x+U.y*a.V.x+U.z*a.W.x, U.x*a.U.y+U.y*a.V.y+U.z*a.W.y, U.x*a.U.z+U.y*a.V.z+U.z*a.W.z,
			V.x*a.U.x+V.y*a.V.x+V.z*a.W.x, V.x*a.U.y+V.y*a.V.y+V.z*a.W.y, V.x*a.U.z+V.y*a.V.z+V.z*a.W.z,
			W.x*a.U.x+W.y*a.V.x+W.z*a.W.x, W.x*a.U.y+W.y*a.V.y+W.z*a.W.y, W.x*a.U.z+W.y*a.V.z+W.z*a.W.z ); }
	// Determinant
	inline _T determinant() const
		{ return U.x*V.y*W.z + U.y*V.z*W.x + U.z*V.x*W.y
				-U.z*V.y*W.x - U.x*V.z*W.y - U.y*V.x*W.z; }
	// Transpose
	inline mat33<_T> transpose() const
		{ return mat33<_T>(
			U.x, V.x, W.x,
			U.y, V.y, W.y,
			U.z, V.z, W.z ); }
	//Trace
	inline _T trace() const
	{ return U.x+V.y+W.z;
		}

	//divergence ,散度（▽·），20161107,
	inline _T divergence() const
	{ return 
	U.x*U.y+V.x*V.z+W.y*W.z;
	}




	// << and >>
	friend inline std::ostream& operator<<( std::ostream& os, const mat33<_T>& a )
		{ os << a.U << std::endl << a.V << std::endl << a.W << std::endl; return os; }
	friend inline std::istream& operator>>( std::istream& is, mat33<_T>& a )
		{ is >> a.U >> a.V >> a.W; return is; }
	friend inline std::ofstream& operator<<( std::ofstream& os, const mat33<_T>& a )
		{ os << a.U << std::endl << a.V << std::endl << a.W << std::endl; return os; }
	friend inline std::ifstream& operator>>( std::ifstream& is, mat33<_T>& a )
		{ is >> a.U >> a.V >> a.W; return is; }

	// Pointer
	inline _T* ptr()
	{ return &U.x; }
	inline const _T* ptr() const
	{ return &U.x; }

	// Data
	static const mat33<_T> E; // Identity matrix
	static const mat33<_T> O; // Zero matrix
};

template<typename _T>
const mat33<_T> mat33<_T>::E=mat33<_T>(
	_T(1),_T(0),_T(0),
	_T(0),_T(1),_T(0),
	_T(0),_T(0),_T(1) );
template<typename _T>
const mat33<_T> mat33<_T>::O=mat33<_T>(
	_T(0),_T(0),_T(0),
	_T(0),_T(0),_T(0),
	_T(0),_T(0),_T(0) );

template<typename _T>
inline mat33<_T> vec3<_T>::pro_transpose(const vec3<_T>& a) const
{ return  mat33<_T>(
	x*a.x, x*a.y, x*a.z,
	y*a.x, y*a.y, y*a.z,
	z*a.x, z*a.y, z*a.z ); }
template<typename _T>
inline mat33<_T> vec3<_T>::pro_trans_self() const
{ return  mat33<_T>(
	x*x, x*y, x*z,
	y*x, y*y, y*z,
	z*x, z*y, z*z ); }
template<typename _T>
inline vec3<_T> vec3<_T>::pro_matrix_left(const mat33<_T>& a) const
{ return  vec3<_T>(
	a.U.x*x+a.U.y*y+a.U.z*z,
	a.V.x*x+a.V.y*y+a.V.z*z,
	a.W.x*x+a.W.y*y+a.W.z*z ); }
template<typename _T>
inline vec3<_T> vec3<_T>::pro_matrix_right(const mat33<_T>& a) const
{ return  vec3<_T>(
	x*a.U.x+y*a.V.x+z*a.W.x,
	x*a.U.y+y*a.V.y+z*a.W.y,
	x*a.U.z+y*a.V.z+z*a.W.z ); }



#endif // #ifndef VEC23_H_


