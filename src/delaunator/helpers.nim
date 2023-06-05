# TODO: inline things
#       iterators sans ids?

import std/[math, sequtils, sets, sugar, options]

import ../delaunator
import clip

# Region clipping based on https://observablehq.com/@mbostock/to-infinity-and-back-again

#[
bounds
- halfedgeIdsOfTriange
- triangleIdOfEdge
- nextHalfedge
- prevHalfedge
- pointIdsOfTriangle
- triangleIdsAdjacentToTriangle
- iterPoints
- iterHullPoints
- iterHullEdges
- iterTriangleEdges
- iterTriangles
- triangleCircumenter
- iterVoronoiEdges
- edgeIdsAroundPoint
- iterVoronoiRegions
- voronoiRegion

etc...
onion
centroids (of regions)
spanning tree
neighbor sites
circles (largest circle fitting the region of each site centered at site)
nearest site
]#


func halfedgeIdsOfTriangle*(t: int32): array[3, int32] =
  ## The halfedge ids of a triangle with id `t`.
  return [3 * t, 3 * t + 1, 3 * t + 2]


func triangleIdOfEdge*(e: int32): int32 =
  ## The id of the triangle for which halfedge with id `e` is a part. (also the id of the 1st point?)
  return floorDiv[int32](e, 3)


func nextHalfedge*(e: int32): int32 =
  ## The id of the next halfedge of the triangle for which halfedge with id `e` is a part.
  if e %% 3 == 2:
    return e - 2
  else:
    return e + 1


func prevHalfedge*(e: int32): int32 =
  # The id of the previous halfedge of the triangle for which halfedge with id `e` is a part.
  if e %% 3 == 0:
    return e + 2
  else:
    return e - 1


iterator iterTriangleEdges*[T](d: Delaunator[T]): tuple[e: int; p, q: array[2, T]] =
  ## Provides an iterator yielding values for each edge of the triangulation.
  ## The values yielded are the id of the halfedge chosen for the edge, an
  ## array describing the point the edge starts at, and an array describing
  ## the point the edge ends at.
  var e = 0
  while e < d.triangles.len:
    if e > d.halfedges[e]:
      let
        pid = d.triangles[e]
        qid = d.triangles[nextHalfedge(int32(e))]
        p = [d.coords[2 * pid], d.coords[2 * pid + 1]]
        q = [d.coords[2 * qid], d.coords[2 * qid + 1]]
      yield (e, p, q)
    inc e


func pointIdsOfTriangle*(d: Delaunator, t: int32): seq[int32] =
  ## The point ids composing the triangle with id `t`.
  return map(halfedgeIdsOfTriangle(t), proc(h: int32): int32 = int32(d.triangles[h]))


iterator iterTriangles*[T](d: Delaunator[T]): tuple[t: int; p1, p2, p3: array[2, T]] =
  ## Provides an iterator yielding values for each triangle of the triangulation.
  ## The values yielded are the id of the triangle, and three arrays, each
  ## describing a point of the triangle.
  var t = 0
  while t < floorDiv(d.triangles.len, 3):
    let
      pids = pointIdsOfTriangle(d, int32(t))
      p1 = [d.coords[2 * pids[0]], d.coords[2 * pids[0] + 1]]
      p2 = [d.coords[2 * pids[1]], d.coords[2 * pids[1] + 1]]
      p3 = [d.coords[2 * pids[2]], d.coords[2 * pids[2] + 1]]
    yield (t, p1, p2, p3)
    inc t


func triangleIdsAdjacentToTriangle*(d: Delaunator, t: int32): seq[int32] =
  ## The triangle ids adjacent to triangle with id `t`. A seq of length 2 or 3.
  var tids = collect(newSeqOfCap(3)):
    for h in halfedgeIdsOfTriangle(t):
      let opposite = d.halfedges[h]
      if opposite >= 0: triangleIdOfEdge(opposite)
  return tids


func triangleCentroid*[T](d: Delaunator[T], t: int32): array[2, T] =
  ## The centroid of triangle with id `t`.
  let
    pids = pointIdsOfTriangle(d, t)
    p1 = [d.coords[2 * pids[0]], d.coords[2 * pids[0] + 1]]
    p2 = [d.coords[2 * pids[1]], d.coords[2 * pids[1] + 1]]
    p3 = [d.coords[2 * pids[2]], d.coords[2 * pids[2] + 1]]
    x = (p1[0] + p2[0] + p3[0]) / 3.0
    y = (p1[1] + p2[1] + p3[1]) / 3.0
  return [x, y]


# TODO: polygonCentroid from clip


func triangleCircumcenter*[T](d: Delaunator[T], t: int32): array[2, T] =
  ## The circumcenter of triangle with id `t`.
  let
    pids = pointIdsOfTriangle(d, t)
    p1 = [d.coords[2 * pids[0]], d.coords[2 * pids[0] + 1]]
    p2 = [d.coords[2 * pids[1]], d.coords[2 * pids[1] + 1]]
    p3 = [d.coords[2 * pids[2]], d.coords[2 * pids[2] + 1]]
    (cx, cy) = circumcenter(p1[0], p1[1], p2[0], p2[1], p3[0], p3[1])
  return [cx, cy]


func edgeIdsAroundPoint*(d: Delaunator, e: int32): seq[int32] =
  ## The ids of all halfedges pointing to the point that halfedge `e` points to.
  var
    start = e
    incomming = start
    edgeIds = newSeqOfCap[int32](10) # TODO: Not sure what's a good starting capacity
  while true:
    edgeIds.add(incomming)
    incomming = d.halfedges[nextHalfedge(incomming)]
    if incomming == -1 or incomming == start: break
  return edgeIds


iterator iterVoronoiEdges*[T](d: Delaunator[T]): tuple[e: int; p, q: array[2, T]] =
  ## Provides an iterator yielding values for each bisecting voronoi edge of the
  ## graph dual to the triangulation. For finite edges, the values yielded are
  ## the id of the halfedge chosen for the bisected edge, an array describing
  ## the circumcenter of the triangle for which that halfedge is a part, and an
  ## array describing the circumcenter of the adjacent triangle for which that
  ## halfedge's compliment is a part. For infinite edges, the ids assigned are
  ## arbitrarily negative. The first array describes the circumcenter of the
  ## hull triangle, and the second the point projected by the ray onto the
  ## delaunator object's defined bounds.
  # First yield edges of finite regions
  var e = 0
  while e < d.triangles.len:
    if e < d.halfedges[e]: # excludes halfedges on hull
      let
        p = triangleCircumcenter[T](d, triangleIdOfEdge(int32(e)))
        q = triangleCircumcenter[T](d, triangleIdOfEdge(d.halfedges[e]))
      yield (e, p, q)
    inc e
  var i = -1 # projected edge ids have negative id, tradeoff :/
  # Second, yield an edge representing one projected
  # ray of each infinite region on the hull.
  for p in d.hull:
    let
      incomming = pointToLeftmostHalfedge(d, int32(p))
      edgeIds = edgeIdsAroundPoint(d, incomming)
      triangleIds = map(edgeIds, proc(h: int32): int32 = triangleIdOfEdge(h))
      vertices = map(triangleIds, proc(t: int32): array[2, T] = triangleCircumcenter(d, t))
      v = p * 4
      pjctd = project(vertices[0], [d.vectors[v], d.vectors[v + 1]], d.bounds.minX, d.bounds.minY, d.bounds.maxX, d.bounds.maxY)
    if pjctd.isSome:
      yield (i, vertices[0], pjctd.get)
      dec i


# TODO: return tuple p as int?
# TODO: // degenerate case? (1 valid point: return the box)
#    if (i === 0 && this.delaunay.hull.length === 1) {
#      return [this.xmax, this.ymin, this.xmax, this.ymax, this.xmin, this.ymax, this.xmin, this.ymin];
#    }
iterator iterVoronoiRegions*[T](d: Delaunator[T]): tuple[p: uint32, verts: seq[array[2, T]]] =
  ## Provides an iterator yielding values for each region of the voronoi diagram.
  ## The values yielded are the id of the point to which the region belongs, and
  ## a seq of vertices describing the region's polygon. Infinite regions are
  ## clipped to the bounds as set on the delaunator object, defaulting to the
  ## minimum and maximum extents of the points provided. Clipped regions which
  ## are completely out of bounds are not yielded, wheras finite regions are
  ## always yielded whether in bounds or not, and are not clipped.
  var
    seen: HashSet[uint32]
    e = 0
  while e < d.triangles.len:
    let p = d.triangles[nextHalfedge(int32(e))]
    if not seen.contains(p):
      seen.incl(p)
      let
        incomming = pointToLeftmostHalfedge(d, int32(p))
        edgeIds = edgeIdsAroundPoint(d, incomming)
        triangleIds = map(edgeIds, proc(h: int32): int32 = triangleIdOfEdge(h))
        vertices = map(triangleIds, proc(t: int32): array[2, T] = triangleCircumcenter(d, t))
        inHullAt = find(d.hull, p)
      if inHullAt == -1: # Not a hull site
        yield (p, vertices)
      else:
        # must clip infinite region
        let
          v = p * 4 # index into ray vectors
          ply = InfConvexPoly[T](points: vertices, v0: [d.vectors[v], d.vectors[v + 1]], vn: [d.vectors[v + 2], d.vectors[v + 3]])
          clipped = clipInfinite[T](ply, d.bounds.minX, d.bounds.minY, d.bounds.maxX, d.bounds.maxY)
        if clipped.len != 0:
          yield (p, clipped)
    inc e


# TODO: return tuple p as int?
# TODO: // degenerate case (1 valid point: return the box)
#    if (i === 0 && this.delaunay.hull.length === 1) {
#      return [this.xmax, this.ymin, this.xmax, this.ymax, this.xmin, this.ymax, this.xmin, this.ymin];
#    }
# TODO: doc that this will return empty verts when completely clipped out of bounds.
proc voronoiRegion*[T](d: Delaunator[T], p: int32): tuple[p: uint32, verts: seq[array[2, T]]] =
  let
    incomming = pointToLeftmostHalfedge(d, p)
    edgeIds = edgeIdsAroundPoint(d, incomming)
    triangleIds = map(edgeIds, proc(h: int32): int32 = triangleIdOfEdge(h))
    vertices = map(triangleIds, proc(t: int32): array[2, T] = triangleCircumcenter(d, t))
    inHullAt = find(d.hull, uint32(p))
  if inHullAt == -1: # Not a hull site
    return (uint32(p), vertices)
  else:
    # must clip infinite region
    let
      v = p * 4 # index into ray vectors
      ply = InfConvexPoly[T](points: vertices, v0: [d.vectors[v], d.vectors[v + 1]], vn: [d.vectors[v + 2], d.vectors[v + 3]])
    return (uint32(p), clipInfinite[T](ply, d.bounds.minX, d.bounds.minY, d.bounds.maxX, d.bounds.maxY))


iterator iterPoints*[T](d: Delaunator[T]): tuple[id: int, p: array[2, T]] =
  ## Provides an iterator yielding values for each point of the triangulation.
  ## The values yielded are the id of the point, and an array describing the
  ## point's location.
  var p = 0
  while p < floorDiv(d.coords.len, 2):
    yield (p, [d.coords[2 * p], d.coords[2 * p + 1]])
    inc p


iterator iterHullPoints*[T](d: Delaunator[T]): tuple[id: int, p: array[2, T]] =
  ## Provides an iterator yielding values for each point of the triangulation's
  ## hull. The values yielded are the id of the point in `d.hull`, and an array
  ## describing the point's location.
  var p = 0
  while p < d.hull.len:
    yield (p, [d.coords[2 * d.hull[p]], d.coords[2 * d.hull[p] + 1]])
    inc p


iterator iterHullEdges*[T](d: Delaunator[T]): tuple[e: int; p, q: array[2, T]] =
  ## Provides an iterator yielding values for each edge of the triangulation's hull.
  ## The values yielded are the id of the halfedge running from p to q, an array
  ## describing the point the edge starts at, and an array describing the point
  ## the edge ends at.
  var e = 0
  while e < d.hull.len:
    let
      pid = d.hull[e]
      qid = hullNext(d, pid)
      p = [d.coords[2 * pid], d.coords[2 * pid + 1]]
      q = [d.coords[2 * qid], d.coords[2 * qid + 1]]
    yield (e, p, q)
    inc e


#TODO: hullPointVectors?
