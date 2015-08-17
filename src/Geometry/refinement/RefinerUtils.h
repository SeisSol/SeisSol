#ifndef _REFINER_UTILS_H_
#define _REFINER_UTILS_H_

#include <cmath>
#include <vector>

namespace refinement {

//------------------------------------------------------------------------------

template<class T>
struct Vec3D {
    T x, y, z;

    Vec3D() : x(0), y(0), z(0) {};
    Vec3D(T X, T Y, T Z) : x(X), y(Y), z(Z) {};
    Vec3D(const T a[3]) : x(a[0]), y(a[1]), z(a[2]) {};

    const T norm() const
    { 
        return std::sqrt(x*x + y*y + z*z);
    }

    Vec3D<T>& operator*=(T c)
    {
        x *= c;
        y *= c;
        z *= c;
        return *this;
    }

    const Vec3D<T> operator-(const Vec3D<T>& b) const
    {
        return Vec3D<T>(x-b.x, y-b.y, z-b.z);
    }

    const Vec3D<T> operator+(const Vec3D<T>& b) const
    {
        return Vec3D<T>(x+b.x, y+b.y, z+b.z);
    }

    const Vec3D<T> operator/(T c) const
    {
        return Vec3D<T>(x/c, y/c, z/c);
    }

    bool operator==(const Vec3D<T> b) const
    {
        return (x==b.x) && (y==b.y) && (z==b.z);
    }

    void dumpData(T*& dest) const
    {
        *dest++ = x;
        *dest++ = y;
        *dest++ = z;
    }
};

//------------------------------------------------------------------------------

template<class T>
const Vec3D<T> middle(const Vec3D<T>& a, const Vec3D<T>& b)
{
    return (a+b)/2;
}

//------------------------------------------------------------------------------

template<class T>
struct Tetrahedron {
    Vec3D<T> a, b, c, d;

    Tetrahedron() {};

    Tetrahedron(const Vec3D<T>& A, const Vec3D<T>& B, const Vec3D<T>& C,
            const Vec3D<T>& D) : a(A), b(B), c(C), d(D) {};

    Tetrahedron(const Vec3D<T>* A, const Vec3D<T>* B, const Vec3D<T>* C,
            const Vec3D<T>* D) : a(*A), b(*B), c(*C), d(*D) {};

    Tetrahedron(const T A[3], const T B[3], const T C[3], const T D[3])
            : a(A), b(B), c(C), d(D) {};


    static const Tetrahedron<T> unitTetrahedron()
    {
        return Tetrahedron(
                Vec3D<T>(0,0,0),
                Vec3D<T>(1,0,0),
                Vec3D<T>(0,1,0),
                Vec3D<T>(0,0,1)
                );
    }

    const Vec3D<T> center() const
    {
        return middle(middle(a, b), middle(c, d));
    }

    void dumpData(T*& dest) const
    {
        a.dumpData(dest);
        b.dumpData(dest);
        c.dumpData(dest);
        d.dumpData(dest);
    }
};

//------------------------------------------------------------------------------

template<class T>
class TetrahedronRefiner {
public:
    virtual void operator()(const Tetrahedron<T>& in, Tetrahedron<T>* out) const = 0;
    virtual unsigned int getDivisionCount() const = 0;
};

//------------------------------------------------------------------------------

template<class T>
class DivideTetrahedronBy4 : public TetrahedronRefiner<T> {
public:
    void operator()(const Tetrahedron<T>& in, Tetrahedron<T>* out) const
    {
        Vec3D<T> mid = in.center();

        *(out++) = Tetrahedron<T>(in.a, in.b, in.c, mid);
        *(out++) = Tetrahedron<T>(in.a, in.b, in.d, mid);
        *(out++) = Tetrahedron<T>(in.a, in.c, in.d, mid);
        *(out++) = Tetrahedron<T>(in.b, in.c, in.d, mid);
    }

    unsigned int getDivisionCount() const
    {
        return 4;
    }
};

//------------------------------------------------------------------------------

template<class T>
class DivideTetrahedronBy8 : public TetrahedronRefiner<T> {
public:
    void operator()(const Tetrahedron<T>& in, Tetrahedron<T>* out) const
    {
        const Vec3D<T>& a = in.a;
        const Vec3D<T>& b = in.b;
        const Vec3D<T>& c = in.c;
        const Vec3D<T>& d = in.d;
        const Vec3D<T> ab = middle(a, b);
        const Vec3D<T> ac = middle(a, c);
        const Vec3D<T> ad = middle(a, d);
        const Vec3D<T> bc = middle(b, c);
        const Vec3D<T> bd = middle(b, d);
        const Vec3D<T> cd = middle(c, d);

        // Corner Cells
        *(out++) = Tetrahedron<T>( a, ab, ac, ad);
        *(out++) = Tetrahedron<T>( b, ab, bc, bd);
        *(out++) = Tetrahedron<T>( c, ac, bc, cd);
        *(out++) = Tetrahedron<T>( d, ad, bd, cd);
        // Inner upper cells
        *(out++) = Tetrahedron<T>(ab, ac, ad, bd);
        *(out++) = Tetrahedron<T>(ab, ac, bc, bd);
        // Inner lower cells
        *(out++) = Tetrahedron<T>(ac, ad, bd, cd);
        *(out++) = Tetrahedron<T>(ac, bc, bd, cd);
    }

    unsigned int getDivisionCount() const
    {
        return 8;
    }
};

//------------------------------------------------------------------------------

template<class T>
class DivideTetrahedronBy32 : public TetrahedronRefiner<T> {
private:
    DivideTetrahedronBy4<T> div4;
    DivideTetrahedronBy8<T> div8;

public:
    void operator()(const Tetrahedron<T>& in, Tetrahedron<T>* out) const
    {
        std::vector<Tetrahedron<T> > buffer(8);
        div8(in, buffer.data());
        for (int i = 0; i < 8; i++)
            div4(buffer[i], out+i*4);
    }

    unsigned int getDivisionCount() const
    {
        return 32;
    }
};

//------------------------------------------------------------------------------

} // namespace

#endif // _REFINER_UTILS_H_
