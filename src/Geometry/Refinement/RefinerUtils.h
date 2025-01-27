// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_GEOMETRY_REFINEMENT_REFINERUTILS_H_
#define SEISSOL_SRC_GEOMETRY_REFINEMENT_REFINERUTILS_H_

#include <vector>

#include <Eigen/Dense>

namespace seissol {
namespace refinement {

//------------------------------------------------------------------------------

template <class T>
const Eigen::Matrix<T, 3, 1> middle(const Eigen::Matrix<T, 3, 1>& a,
                                    const Eigen::Matrix<T, 3, 1>& b) {
  return (a + b) / static_cast<T>(2);
}

//------------------------------------------------------------------------------

template <class T>
struct Tetrahedron {
  /** Indices of the vertices */
  unsigned int i, j, k, l;
  Eigen::Matrix<T, 3, 1> a, b, c, d;

  Tetrahedron() : i(0), j(0), k(0), l(0) {};

  Tetrahedron(const Eigen::Matrix<T, 3, 1>& A,
              const Eigen::Matrix<T, 3, 1>& B,
              const Eigen::Matrix<T, 3, 1>& C,
              const Eigen::Matrix<T, 3, 1>& D)
      : i(0), j(0), k(0), l(0), a(A), b(B), c(C), d(D) {};

  Tetrahedron(const Eigen::Matrix<T, 3, 1>& A,
              const Eigen::Matrix<T, 3, 1>& B,
              const Eigen::Matrix<T, 3, 1>& C,
              const Eigen::Matrix<T, 3, 1>& D,
              unsigned int I,
              unsigned int J,
              unsigned int K,
              unsigned int L)
      : i(I), j(J), k(K), l(L), a(A), b(B), c(C), d(D) {};

  Tetrahedron(const T A[3],
              const T B[3],
              const T C[3],
              const T D[3],
              unsigned int I,
              unsigned int J,
              unsigned int K,
              unsigned int L)
      : i(I), j(J), k(K), l(L), a(Eigen::Matrix<T, 3, 1>(A)), b(Eigen::Matrix<T, 3, 1>(B)),
        c(Eigen::Matrix<T, 3, 1>(C)), d(Eigen::Matrix<T, 3, 1>(D)) {};

  static const Tetrahedron<T> unitTetrahedron() {
    return Tetrahedron(Eigen::Matrix<T, 3, 1>(0, 0, 0),
                       Eigen::Matrix<T, 3, 1>(0, 0, 1),
                       Eigen::Matrix<T, 3, 1>(1, 0, 0),
                       Eigen::Matrix<T, 3, 1>(0, 1, 0));
  }

  const Eigen::Matrix<T, 3, 1> center() const { return middle(middle(a, b), middle(c, d)); }
};

//------------------------------------------------------------------------------

template <class T>
class TetrahedronRefiner {
  public:
  virtual ~TetrahedronRefiner() {}
  /** Generates the new cells */
  virtual void refine(const Tetrahedron<T>& in,
                      unsigned int addVertexStart,
                      Tetrahedron<T>* out,
                      Eigen::Matrix<T, 3, 1>* addVertices) const = 0;
  /** The additional number of vertices generated per original tetrahedron */
  virtual unsigned int additionalVerticesPerCell() const = 0;
  virtual unsigned int getDivisionCount() const = 0;
};

//------------------------------------------------------------------------------

template <class T>
class IdentityRefiner : public TetrahedronRefiner<T> {
  public:
  void refine(const Tetrahedron<T>& in,
              unsigned int addVertexStart,
              Tetrahedron<T>* out,
              Eigen::Matrix<T, 3, 1>* addVertices) const {
    out[0] = in;
  }

  unsigned int additionalVerticesPerCell() const { return 0; }

  unsigned int getDivisionCount() const { return 1; }
};

//------------------------------------------------------------------------------

template <class T>
class DivideTetrahedronBy4 : public TetrahedronRefiner<T> {
  public:
  void refine(const Tetrahedron<T>& in,
              unsigned int addVertexStart,
              Tetrahedron<T>* out,
              Eigen::Matrix<T, 3, 1>* addVertices) const {
    addVertices[0] = in.center();

    out[0] = Tetrahedron<T>(in.a, in.b, in.c, addVertices[0], in.i, in.j, in.k, addVertexStart);
    out[1] = Tetrahedron<T>(in.a, in.b, in.d, addVertices[0], in.i, in.j, in.l, addVertexStart);
    out[2] = Tetrahedron<T>(in.a, in.c, in.d, addVertices[0], in.i, in.k, in.l, addVertexStart);
    out[3] = Tetrahedron<T>(in.b, in.c, in.d, addVertices[0], in.j, in.k, in.l, addVertexStart);
  }

  unsigned int additionalVerticesPerCell() const { return 1; }

  unsigned int getDivisionCount() const { return 4; }
};

//------------------------------------------------------------------------------

template <class T>
class DivideTetrahedronBy8 : public TetrahedronRefiner<T> {
  public:
  void refine(const Tetrahedron<T>& in,
              unsigned int addVertexStart,
              Tetrahedron<T>* out,
              Eigen::Matrix<T, 3, 1>* addVertices) const {
    const Eigen::Matrix<T, 3, 1>& a = in.a;
    const Eigen::Matrix<T, 3, 1>& b = in.b;
    const Eigen::Matrix<T, 3, 1>& c = in.c;
    const Eigen::Matrix<T, 3, 1>& d = in.d;
    addVertices[0] = middle(a, b);
    addVertices[1] = middle(a, c);
    addVertices[2] = middle(a, d);
    addVertices[3] = middle(b, c);
    addVertices[4] = middle(b, d);
    addVertices[5] = middle(c, d);

    const Eigen::Matrix<T, 3, 1>& ab = addVertices[0];
    const Eigen::Matrix<T, 3, 1>& ac = addVertices[1];
    const Eigen::Matrix<T, 3, 1>& ad = addVertices[2];
    const Eigen::Matrix<T, 3, 1>& bc = addVertices[3];
    const Eigen::Matrix<T, 3, 1>& bd = addVertices[4];
    const Eigen::Matrix<T, 3, 1>& cd = addVertices[5];

    const unsigned int iab = addVertexStart;
    const unsigned int iac = addVertexStart + 1;
    const unsigned int iad = addVertexStart + 2;
    const unsigned int ibc = addVertexStart + 3;
    const unsigned int ibd = addVertexStart + 4;
    const unsigned int icd = addVertexStart + 5;

    out[0] = Tetrahedron<T>(a, ab, ac, ad, in.i, iab, iac, iad);
    out[1] = Tetrahedron<T>(b, ab, bc, bd, in.j, iab, ibc, ibd);
    out[2] = Tetrahedron<T>(c, ac, bc, cd, in.k, iac, ibc, icd);
    out[3] = Tetrahedron<T>(d, ad, bd, cd, in.l, iad, ibd, icd);
    // Inner upper cells
    out[4] = Tetrahedron<T>(ab, ac, ad, bd, iab, iac, iad, ibd);
    out[5] = Tetrahedron<T>(ab, ac, bc, bd, iab, iac, ibc, ibd);
    // Inner lower cells
    out[6] = Tetrahedron<T>(ac, ad, bd, cd, iac, iad, ibd, icd);
    out[7] = Tetrahedron<T>(ac, bc, bd, cd, iac, ibc, ibd, icd);
  }

  unsigned int additionalVerticesPerCell() const { return 6; }

  unsigned int getDivisionCount() const { return 8; }
};

//------------------------------------------------------------------------------

template <class T>
class DivideTetrahedronBy32 : public TetrahedronRefiner<T> {
  private:
  DivideTetrahedronBy4<T> div4;
  DivideTetrahedronBy8<T> div8;

  public:
  void refine(const Tetrahedron<T>& in,
              unsigned int addVertexStart,
              Tetrahedron<T>* out,
              Eigen::Matrix<T, 3, 1>* addVertices) const {
    auto tmp = std::vector<Tetrahedron<T>>(div8.getDivisionCount());

    div8.refine(in, addVertexStart, tmp.data(), addVertices);

    addVertexStart += div8.additionalVerticesPerCell();
    addVertices += div8.additionalVerticesPerCell();

    for (unsigned int i = 0; i < div8.getDivisionCount(); i++) {
      div4.refine(tmp[i],
                  addVertexStart + i * div4.additionalVerticesPerCell(),
                  out + (i * div4.getDivisionCount()),
                  addVertices + (i * div4.additionalVerticesPerCell()));
    }
  }

  unsigned int additionalVerticesPerCell() const {
    return div8.additionalVerticesPerCell() +
           div8.getDivisionCount() * div4.additionalVerticesPerCell();
  }

  unsigned int getDivisionCount() const {
    return div4.getDivisionCount() * div8.getDivisionCount();
  }
};

//------------------------------------------------------------------------------

} // namespace refinement
} // namespace seissol

#endif // SEISSOL_SRC_GEOMETRY_REFINEMENT_REFINERUTILS_H_
