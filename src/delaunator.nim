
import std/math
from std/fenv import epsilon


var EDGE_STACK: array[512, uint32]

type
  Delaunator*[T] = ref object
    coords*: seq[T]
    triangles*: seq[uint32]
    halfedges*: seq[int32]
    hull*: seq[uint32]

    # Arrays that will store the triangulation graph
    trianglesLen: int
    d_triangles: seq[uint32]
    d_halfedges: seq[int32]

    # Temporary arrays for tracking the edges of the advancing convex hull
    d_hashSize: int
    d_hullPrev: seq[uint32]  # edge to prev edge
    d_hullNext: seq[uint32]  # edge to next edge
    d_hullTri:  seq[uint32]  # edge to adjacent triangle
    d_hullHash: seq[int32]   # angular edge hash

    # Temporary arrays for sorting points
    d_ids:   seq[uint32]
    d_dists: seq[T]


func defaultGetX[P, T](p: P): T = p[0]
func defaultGetY[P, T](p: P): T = p[1]


proc update(this: var Delaunator) =
  this.triangles = this.d_triangles[1 .. ^1]
  discard



proc fromCoords*[T](coordinates: seq[T]): Delaunator[T] =
  result = Delaunator[T](coords: coordinates)
  let
    n = ashr(coordinates.len, 1) # n points
    maxTriangles = max(2 * n - 5, 0)
  result.trianglesLen = 0
  result.d_triangles = newSeq[uint32](maxTriangles * 3)
  result.d_halfedges = newSeq[int32](maxTriangles * 3)
  result.d_hashSize = ceil(sqrt(n.toFloat)).toInt
  result.d_hullPrev = newSeq[uint32](n)
  result.d_hullNext = newSeq[uint32](n)
  result.d_hullTri = newSeq[uint32](n)
  result.d_hullHash = newSeq[int32](result.d_hashSize)
  result.d_ids = newSeq[uint32](n)
  result.d_dists = newSeq[T](n)
  update(result)


proc fromPoints*[P, T](points: seq[P]; getX: proc (p: P): T = defaultGetX; getY: proc (p: P): T = defaultGetY): Delaunator[T] =
  var
    coords = newSeq[T](points.len * 2)
  for i, point in points:
      coords[2 * i] = getX(point)
      coords[2 * i + 1] = getY(point)
  fromCoords(coords)


when isMainModule:
  type
    Point = tuple[x, y: float64]

  var
    d: Delaunator[float64]
    points: seq[Point]

  #points = @[(1.0'f32, 1.0'f32), (-1.0'f32, 1.0'f32), (0.0'f32, -1.0'f32)]
  points = @[(1.0, 1.0), (-1.0, 1.0), (0.0, -1.0)]
  d = fromPoints[Point, float64](points)

  echo repr(d)
  echo "Triangles: ", d.triangles
  echo "Halfedges: ", d.halfedges
