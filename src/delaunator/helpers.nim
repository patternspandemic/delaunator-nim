# TODO: inline things
#       iterators sans ids?

import std/[sequtils, sets, sugar]


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

midpoint
perpBisectSlope
yIntercept
abcOfBisector

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


iterator iterTriangleEdges*(d: Delaunator): tuple[e: int32; p, q: array[2, SomeFloat]] =
  ## Provides an iterator yielding values for each edge of the triangulation.
  ## The values yielded are the id of the halfedge chosen for the edge, an
  ## array describing the point the edge starts at, and an array describing
  ## the point the edge ends at.
  var e = 0
  while true:
    if e > d.halfedges[e]:
      let
        pid = d.triangles[e]
        qid = d.triangles[nextHalfedge(e)]
        p: array[2, SomeFloat] = d.coords[(2 * pid) ..< (2 * pid + 2)]
        q: array[2, SomeFloat] = d.coords[(2 * qid) ..< (2 * qid + 2)]
      yield (e, p, q)
    if not e < d.triangles.len: break
    inc e


func pointIdsOfTriangle*(d: Delaunator, t: int32): seq[int32] =
  ## The point ids composing the triangle with id `t`.
  return map(halfedgeIdsOfTriangle(t), proc(h: int32): int32 = d.triangles[h])


iterator iterTriangles*(d: Delaunator): tuple[t: int32; p1, p2, p3: array[2, SomeFloat]] =
  ## Provides an iterator yielding values for each triangle of the triangulation.
  ## The values yielded are the id of the triangle, and three arrays, each
  ## describing a point of the triangle.
  var t = 0
  while true:
    let
      pids = pointIdsOfTriangle(d, t)
      p1: array[2, SomeFloat] = d.coords[(2 * pids[0]) ..< (2 * pids[0] + 2)]
      p2: array[2, SomeFloat] = d.coords[(2 * pids[1]) ..< (2 * pids[1] + 2)]
      p3: array[2, SomeFloat] = d.coords[(2 * pids[2]) ..< (2 * pids[2] + 2)]
    yield (t, p1, p2, p3)
    if not t < (d.triangles.len / 3): break
    inc t


func triangleIdsAdjacentToTriangle*(d: Delaunator, t: int32): seq[int32] =
  ## The triangle ids adjacent to triangle with id `t`. A seq of length 2 or 3.
  var tids = collect(newSeqOfCap(3)):
    for h in halfedgeIdsOfTriangle(t):
      let opposite = d.halfedges[h]
      if opposite >= 0: triangleIdOfEdge(opposite)
  return tids


func triangleCircumcenter*(d: Delaunator, t: int32): array[2, SomeFLoat] =
  ## The circumcenter of triangle with id `t`.
  let
    pids = pointIdsOfTriangle(d, t)
    p1: array[2, SomeFloat] = d.coords[(2 * pids[0]) ..< (2 * pids[0] + 2)]
    p2: array[2, SomeFloat] = d.coords[(2 * pids[1]) ..< (2 * pids[1] + 2)]
    p3: array[2, SomeFloat] = d.coords[(2 * pids[2]) ..< (2 * pids[2] + 2)]
    (cx, cy) = circumcenter(p1[0], p1[1], p2[0], p2[1], p3[0], p3[1])
  return [cx, cy]

iterator iterVoronoiEdges*(d: Delaunator): tuple[e: int32; p, q: array[2, SomeFloat]] =
  ## Provides an iterator yielding values for each bisecting voronoi edge of the
  ## graph dual to the triangulation. The values yielded are the id of the
  ## halfedge chosen for the bisected edge, an array describing the
  ## circumcenter of the triangle for which that halfedge is a part, and an
  ## array describing the circumcenter of the adjacent triangle for which that
  ## halfedge's compliment is a part. Excluded by default are infinite voronoi
  ## edges, which result from bisecting halfedges on the hull. For such bisectors
  ## to be included, a bounding region must be provided to clip these edges against.
  ## TODO: Support voronoi edges of the hull. (via overloaded iterator?)
  var e = 0
  while true:
    if e < d.halfedges[e]: # excludes halfedges on hull
      let
        p = triangleCircumcenter(d, triangleIdOfEdge(e))
        q = triangleCircumcenter(d, triangleIdOfEdge(d.halfedges[e]))
      yield (e, p, q)
    if not e < d.triangles.len: break
    inc e


func edgeIdsAroundPoint*(d: Delaunator, e: int32): seq[int32] =
  ## The ids of all halfedges pointing to the point that halfedge `e` points to.
  ## TODO: test whether this works for case when `e` points to a point on the
  ## hull. Should now be okay with setting `start` to leftmost halfedge.
  var
    # Get the leftmost edge pointing to the point e is pointing to.
    start = d.d_pointToLeftmostHalfedgeIndex[nextHalfedge(e)]
    incomming = start
    edgeIds = newSeqOfCap[int32](10) # TODO: Not sure what's a good starting capacity
  while true:
    edgeIds.add(incomming)
    incomming = d.halfedges[nextHalfedge(incomming)]
    if incomming == -1 or incomming == start: break
  return edgeIds


func iterVoronoiRegions*(d: Delaunator): tuple[p: uint32, verts: seq[array[2, SomeFLoat]]] =
  ## Provides an iterator yielding values for each region of the voronoi diagram.
  ## The values yielded are the id of the point to which the region belongs, and
  ## a seq of vertices describing the region's polygon. Excluded by default
  ## are infinite unbound regions. For such regions to be included, a bounding
  ## region must be provided to clip these regions against.
  ## TODO: Support voronoi regions of the hull. (via overloaded iterator?)
  var
    # Do not yield regions of hull points by default
    seen = toHashSet(d.hull)
    e = 0
  while true:
    let p = d.triangles[nextHalfedge(e)]
    if not seen.contains(p):
      seen.incl(p)
      let
        edgeIds = edgeIdsAroundPoint(d, e)
        triangleIds = map(edgeIds, proc(h: int32): int32 = triangleIdOfEdge(h))
        vertices = map(triangleIds, proc(t: int32): array[2, SomeFloat] = triangleCircumcenter(d, t))
      yield (p, vertices)
    if not e < d.triangles.len: break
    inc e


# FIXME: Does not provide clipped regions for points on hull
#        Param p and shadowed p are the same?
func voronoiRegion*(d: Delaunator, p: uint32): tuple[p: uint32, verts: seq[array[2, SomeFloat]]] =
  let
    incomming = d.d_pointToLeftmostHalfedgeIndex[p]
    p = d.triangles[nextHalfedge(incomming)] # FIXME: prob remove
    edgeIds = edgeIdsAroundPoint(d, incomming)
    triangleIds = map(edgeIds, proc(h: int32): int32 = triangleIdOfEdge(h))
    vertices = map(triangleIds, proc(t: int32): array[2, SomeFloat] = triangleCircumcenter(d, t))
  return (p, vertices)


func iterPoints*(d: Delaunator): tuple[id: int32, p: array[2, SomeFloat]] =
  ## Provides an iterator yielding values for each point of the triangulation.
  ## The values yielded are the id of the point, and an array describing the
  ## point's location.
  var p = 0
  while true:
    yield (p, d.coords[(2 * p) ..< (2 * p + 2)])
    if not p < (d.coords.len / 2): break
    inc p


func iterHullPoints*(d: Delaunator): tuple[id: int32, p: array[2, SomeFloat]] =
  ## Provides an iterator yielding values for each point of the triangulation's
  ## hull. The values yielded are the id of the point in `d.hull`, and an array
  ## describing the point's location.
  var p = 0
  while true:
    yield (p, d.coords[(2 * d.hull[p]) ..< (2 * d.hull[p] + 2)])
    if not p < d.hull.len: break
    inc p


# FIXME: Not sure this works
func iterHullEdges*(d: Delaunator): tuple[e: int32; p, q: array[2, SomeFloat]] =
  ## Provides an iterator yielding values for each edge of the triangulation's hull.
  ## The values yielded are the id of the halfedge running from p to q, an array
  ## describing the point the edge starts at, and an array describing the point
  ## the edge ends at.
  var e = 0
  while true:
    let
      pid = d.hull[e]
      qid = d.d_hullNext[pid] # d_hullNext may not be accessable?
      p: array[2, SomeFloat] = d.coords[(2 * pid) ..< (2 * pid + 2)]
      q: array[2, SomeFloat] = d.coords[(2 * qid) ..< (2 * qid + 2)]
    yield (e, p, q)
    if not e < d.hull.len: break
    inc e
