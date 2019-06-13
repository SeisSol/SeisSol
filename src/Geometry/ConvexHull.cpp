#include "ConvexHull.h"
#include "MeshDefinition.h"
#include "MeshTools.h"

#include <algorithm>
#include <cassert>
#include <unordered_set>

#include <iostream> // TODO(Lukas) Remove

bool PointPlane::operator< (const PointPlane &other) const {
  return coords[0] < other.coords[0] ||
	(coords[0] == other.coords[0] && coords[1] < other.coords[1]);
}

bool PointPlane::operator== (const PointPlane &other) const {
  const auto diffX = coords[0] - other.coords[0];
  const auto diffY = coords[1] - other.coords[1];
  const auto diff = diffX * diffX + diffY * diffY;
  const auto eps = 10e-9;
  return diff < eps;
}

double turn(const PointPlane &O,
	    const PointPlane &A,
	    const PointPlane &B) {
  return (A.coords[0] - O.coords[0]) * (B.coords[1] - O.coords[1])
    - (A.coords[1] - O.coords[1]) * (B.coords[0] - O.coords[0]); 
}

// TODO(Lukas) Maybe don't compare with 0 but rather with eps to allow very
// small changes
std::vector<PointPlane> computeConvexHull2D(std::vector<PointPlane> points) {
  auto hullIdx = 0;
  if (points.size() <= 3) {
    // Points already forms a convex hull!
    return points; // TODO(Lukas) Is this the convex hull?
  }

  auto convexHull = std::vector<PointPlane>(2 * points.size());

  // TODO(Lukas): Remove this dirty hack.
  for (auto &point : points) {
    break;
    // Round all points s.t. convex hull looks 'nicer'
    const auto ndigits = 5;
    const auto multp = std::pow(10, ndigits);
    const long long t1 = std::round(point.coords[0] * multp);
    const long long t2 = std::round(point.coords[1] * multp);
    point.coords[0] = t1 / static_cast<double>(multp);
    point.coords[1] = t2 / static_cast<double>(multp);
  }
  std::sort(points.begin(), points.end());

  // Remove duplicate points
  points.erase(std::unique(points.begin(), points.end()), points.end());

  /* lower hull */
  for (int i = 0; i < points.size(); ++i) {
    while (hullIdx >= 2
	   && turn(convexHull[hullIdx-2], convexHull[hullIdx-1], points[i]) <= 0) {
      --hullIdx;
    }
    convexHull[hullIdx++] = points[i];
  }
 
  /* upper hull */
  const auto t = hullIdx + 1;
  for (int i = points.size()-2; i >= 0; --i) {
    while (hullIdx >= t
	   && turn(convexHull[hullIdx-2], convexHull[hullIdx-1], points[i]) <= 0) {
      --hullIdx;
    }
    convexHull[hullIdx++] = points[i];
  }
  

  convexHull.resize(hullIdx-1);

  return convexHull;
}

// TODO(Lukas): Might not work if point is exactly on the convex hull
bool isInsideConvexPolygon(const PointPlane &point,
			   const std::vector<PointPlane> &vertices) {
  const auto eps = 10e-5; // TODO(Lukas) Does eps make sense here?
  
  // Check if all edges make a clockwise turn with our point.
  for (int i = 0; i < vertices.size(); ++i) {
    const auto &coordsA = vertices[i];
    const auto &coordsB = vertices[(i+1) % vertices.size()];
    const auto cross = turn(coordsA, coordsB, point);
    if (cross < eps) return false;
  }
  return true;
}
