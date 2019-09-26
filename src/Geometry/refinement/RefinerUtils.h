/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2015, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 */

#ifndef _REFINER_UTILS_H_
#define _REFINER_UTILS_H_

#include <vector>

#include <glm/vec3.hpp>
#include <glm/gtc/type_ptr.hpp>

namespace seissol {
namespace refinement {

//------------------------------------------------------------------------------

template<class T>
const glm::tvec3<T> middle(const glm::tvec3<T>& a, const glm::tvec3<T>& b)
{
    return (a+b)/static_cast<T>(2);
}

//------------------------------------------------------------------------------

template<class T>
struct Tetrahedron {
	/** Indices of the vertices */
	unsigned int i, j, k, l;
    glm::tvec3<T> a, b, c, d;

    Tetrahedron()
    	: i(0), j(0), k(0), l(0)
	{};

    Tetrahedron(const glm::tvec3<T>& A, const glm::tvec3<T>& B, const glm::tvec3<T>& C,
            const glm::tvec3<T>& D)
    	: i(0), j(0), k(0), l(0),
		  a(A), b(B), c(C), d(D)
    {};

    Tetrahedron(const glm::tvec3<T>& A, const glm::tvec3<T>& B, const glm::tvec3<T>& C,
            const glm::tvec3<T>& D,
			unsigned int I, unsigned int J, unsigned int K, unsigned int L)
    	: i(I), j(J), k(K), l(L),
		  a(A), b(B), c(C), d(D)
    {};

    Tetrahedron(const T A[3], const T B[3], const T C[3], const T D[3],
    		unsigned int I, unsigned int J, unsigned int K, unsigned int L)
            : i(I), j(J), k(K), l(L),
			  a(glm::make_vec3(A)), b(glm::make_vec3(B)),
			  c(glm::make_vec3(C)), d(glm::make_vec3(D))
    {};

    static const Tetrahedron<T> unitTetrahedron()
    {
        return Tetrahedron(
                glm::tvec3<T>(0,0,0),
                glm::tvec3<T>(0,0,1),
                glm::tvec3<T>(1,0,0),
                glm::tvec3<T>(0,1,0));
    }

    const glm::tvec3<T> center() const
    {
        return middle(middle(a, b), middle(c, d));
    }
};

//------------------------------------------------------------------------------

template<class T>
class TetrahedronRefiner {
public:
	virtual ~TetrahedronRefiner() {}
	/** Generates the new cells */
    virtual void refine(const Tetrahedron<T>& in, unsigned int addVertexStart,
    		Tetrahedron<T>* out, glm::tvec3<T>* addVertices) const = 0;
    /** The additional number of vertices generated per original tetrahedron */
    virtual unsigned int additionalVerticesPerCell() const = 0;
    virtual unsigned int getDivisionCount() const = 0;
};

//------------------------------------------------------------------------------

template<class T>
class IdentityRefiner : public TetrahedronRefiner<T> {
public:
	void refine(const Tetrahedron<T>& in, unsigned int addVertexStart,
		Tetrahedron<T>* out, glm::tvec3<T>* addVertices) const
	{
		out[0] = in;
	}

    unsigned int additionalVerticesPerCell() const
    {
    	return 0;
    }

    unsigned int getDivisionCount() const
    {
        return 1;
    }
};

//------------------------------------------------------------------------------

template<class T>
class DivideTetrahedronBy4 : public TetrahedronRefiner<T> {
public:
	void refine(const Tetrahedron<T>& in, unsigned int addVertexStart,
		Tetrahedron<T>* out, glm::tvec3<T>* addVertices) const
	{
		addVertices[0] = in.center();

        out[0] = Tetrahedron<T>(in.a, in.b, in.c, addVertices[0],
        		in.i, in.j, in.k, addVertexStart);
        out[1] = Tetrahedron<T>(in.a, in.b, in.d, addVertices[0],
        		in.i, in.j, in.l, addVertexStart);
        out[2] = Tetrahedron<T>(in.a, in.c, in.d, addVertices[0],
        		in.i, in.k, in.l, addVertexStart);
        out[3] = Tetrahedron<T>(in.b, in.c, in.d, addVertices[0],
        		in.j, in.k, in.l, addVertexStart);
	}

    unsigned int additionalVerticesPerCell() const
    {
    	return 1;
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
	void refine(const Tetrahedron<T>& in, unsigned int addVertexStart,
		Tetrahedron<T>* out, glm::tvec3<T>* addVertices) const
	{
        const glm::tvec3<T>& a = in.a;
        const glm::tvec3<T>& b = in.b;
        const glm::tvec3<T>& c = in.c;
        const glm::tvec3<T>& d = in.d;
        addVertices[0] = middle(a, b);
        addVertices[1] = middle(a, c);
        addVertices[2] = middle(a, d);
        addVertices[3] = middle(b, c);
        addVertices[4] = middle(b, d);
        addVertices[5] = middle(c, d);

        const glm::tvec3<T>& ab = addVertices[0];
        const glm::tvec3<T>& ac = addVertices[1];
        const glm::tvec3<T>& ad = addVertices[2];
        const glm::tvec3<T>& bc = addVertices[3];
        const glm::tvec3<T>& bd = addVertices[4];
        const glm::tvec3<T>& cd = addVertices[5];

        const unsigned int iab = addVertexStart;
        const unsigned int iac = addVertexStart+1;
        const unsigned int iad = addVertexStart+2;
        const unsigned int ibc = addVertexStart+3;
        const unsigned int ibd = addVertexStart+4;
        const unsigned int icd = addVertexStart+5;


        out[0] = Tetrahedron<T>( a, ab, ac, ad,
        		in.i, iab, iac, iad);
        out[1] = Tetrahedron<T>( b, ab, bc, bd,
        		in.j, iab, ibc, ibd);
        out[2] = Tetrahedron<T>( c, ac, bc, cd,
        		in.k, iac, ibc, icd);
        out[3] = Tetrahedron<T>( d, ad, bd, cd,
        		in.l, iad, ibd, icd);
        // Inner upper cells
        out[4] = Tetrahedron<T>(ab, ac, ad, bd,
        		iab, iac, iad, ibd);
        out[5] = Tetrahedron<T>(ab, ac, bc, bd,
        		iab, iac, ibc, ibd);
        // Inner lower cells
        out[6] = Tetrahedron<T>(ac, ad, bd, cd,
        		iac, iad, ibd, icd);
        out[7] = Tetrahedron<T>(ac, bc, bd, cd,
        		iac, ibc, ibd, icd);
	}

    unsigned int additionalVerticesPerCell() const
    {
    	return 6;
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
	void refine(const Tetrahedron<T>& in, unsigned int addVertexStart,
		Tetrahedron<T>* out, glm::tvec3<T>* addVertices) const
	{
	        auto tmp = std::vector<Tetrahedron<T>>(div8.getDivisionCount());

		div8.refine(in, addVertexStart, tmp.data(), addVertices);

		addVertexStart += div8.additionalVerticesPerCell();
		addVertices += div8.additionalVerticesPerCell();

		for (unsigned int i = 0; i < div8.getDivisionCount(); i++) {
			div4.refine(tmp[i],
					addVertexStart + i*div4.additionalVerticesPerCell(),
					out + (i*div4.getDivisionCount()),
					addVertices + (i*div4.additionalVerticesPerCell()));
		}
	}

    unsigned int additionalVerticesPerCell() const
    {
    	return div8.additionalVerticesPerCell()
    			+ div8.getDivisionCount() * div4.additionalVerticesPerCell();
    }

    unsigned int getDivisionCount() const
    {
        return div4.getDivisionCount() * div8.getDivisionCount();
    }
};

//------------------------------------------------------------------------------

} // namespace
}

#endif // _REFINER_UTILS_H_
