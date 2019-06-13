#ifndef CONVEX_HULL_H
#define CONVEX_HULL_H

#include <vector>

// TODO Namespace, code quality, etc.

// Implementation taken from
// https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain#C++

struct PointPlane {
  double coords[2];
  
  bool operator< (const PointPlane &other) const; 

  bool operator== (const PointPlane &other) const;
};

double turn(const PointPlane &O,
	    const PointPlane &A,
	    const PointPlane &B);

std::vector<PointPlane> computeConvexHull2D(std::vector<PointPlane> points);

bool isInsideConvexPolygon(const PointPlane &point,
			   const std::vector<PointPlane> &points);

#endif // CONVEX_HULL_H
