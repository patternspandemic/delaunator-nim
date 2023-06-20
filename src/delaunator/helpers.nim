# Region clipping based on https://observablehq.com/@mbostock/to-infinity-and-back-again

import std/[math, options, sequtils, random, sets, sugar]

import ../delaunator
import clip


#[
etc...
boundsCenter, hullCentroid
onion iterOnionLayers
spanning tree
neighbor sites iterPointNeighbors
circles (largest circle fitting the tri/region centered at some center)
nearest/farthest neighbor
hullPointVectors?
]#

proc normalDist*[F](count: int, mu = 0.0; sigma = 1.0, r: var Rand = randState()): seq[F] =
  ## Returns a `seq` of `count` *normally* distributed random values of given
  ## `mu` and `sigme`.
  result = collect(newSeqOfCap(count)):
    for i in 1 .. count: F(gauss(r, mu, sigma))


proc uniformDist*[F](count: int, a: F, b: F, r: var Rand = randState()): seq[F] =
  ## Returns a `seq` of `count` *uniformly* distributed random values in the
  ## range `a` .. `b`.
  result = collect(newSeqOfCap(count)):
    for i in 1 .. count: rand[F](r, a .. b)


proc grid*[F](
  ## Returns a flattened `seq` of points in an `n` x `m` grid distributed within
  ## a horizontal range `fromN` .. `toN`, and a vertical range `formM` .. `toM`.
  ## By default, points begin and end at range extents, but can be inset by
  ## setting `leftEdge`, `rightEdge`, `topEdge` or `bottomEdge` to false. Set an
  ## amount of horizontal or vertical jitter in the range 0.0 .. 1.0 with `jitterN`
  ## and `jitterM` respectively. When `pinHull` is true, points on the grid's
  ## hull are not jittered.
  n, m: int; # n x m sites
  fromN, toN, fromM, toM: F; # grid range
  leftEdge, rightEdge, topEdge, bottomEdge: bool = true; # sites on edges?
  jitterN, jitterM: F = 0.0; pinHull: bool = false; # jitter
  r: var Rand = randState()
): seq[F] =
  let
    count = n * m
    nWidth = toN - fromN
    mWidth = toM - fromM
  var
    nStepRatio = n - 1
    mStepRatio = m - 1
  if not leftEdge: nStepRatio += 1
  if not rightEdge: nStepRatio += 1
  if not topEdge: mStepRatio += 1
  if not bottomEdge: mStepRatio += 1
  let
    nStep = nWidth / F(nStepRatio)
    mStep = mWidth / F(mStepRatio)
    nStart = fromN + (if not leftEdge: nStep else: F(0.0))
    mStart = fromM + (if not topEdge: mStep else: F(0.0))
  var coords = newSeqOfCap[F](count)
  for i in 0 ..< n:
    for j in 0 ..< m:
      let
        x = nStart + F(i) * nStep
        y = mStart + F(j) * mStep
        xJitter = jitterN * (rand[F](r, nStep) - F(0.5) * nStep)
        yJitter = jitterM * (rand[F](r, mstep) - F(0.5) * mStep)
        pin = if pinHull and ((i == 0 or i == n - 1) or (j == 0 or j == m - 1)): F(0.0)
              else: F(1.0)
      coords.add(x + xJitter * pin)
      coords.add(y + yJitter * pin)
  return coords


proc zipPolarDists*[F](r, t: openArray[F], xOffset, yOffset: F = F(0.0)): seq[F] =
  ## Returns a flattened `seq` of points by zipping two collections of values,
  ## `r` -treated as a radii, and `t` -treated as an angles, converting each pair
  ## of polar coordinates into their cartesian counterpart. Translate the polar
  ## origin with `xOffset` and `yOffset`.
  var
    m = min(r.len, t.len)
    i = 0
  newSeq(result, m * 2)
  while i < m:
    result[i * 2] = r[i] * cos(t[i]) + xOffset
    result[i * 2 + 1] = r[i] * sin(t[i]) + yOffset
    inc i


proc zipDists*[F](s1, s2: openArray[F]): seq[F] =
  ## Returns a flattened `seq` of points by zipping two collections of values,
  ## `s1` and `s2`
  var
    m = min(s1.len, s2.len)
    i = 0
  newSeq(result, m * 2)
  while i < m:
    result[i * 2] = s1[i]
    result[i * 2 + 1] = s2[i]
    inc i


func halfedgeIdsOfTriangle*(tid: uint32): array[3, int32] =
  ## The halfedge ids of a triangle with id `tid`.
  return [3 * int32(tid), 3 * int32(tid) + 1, 3 * int32(tid) + 2]


func triangleIdOfEdge*(eid: int32): uint32 =
  ## The id of the triangle for which halfedge with id `eid` is a part. (also the id of the 1st point?)
  return floorDiv(uint32(eid), 3'u32)


func nextHalfedge*(eid: int32): int32 =
  ## The id of the next halfedge of the triangle for which halfedge with id `eid` is a part.
  if eid %% 3 == 2:
    return eid - 2
  else:
    return eid + 1


func prevHalfedge*(eid: int32): int32 =
  ## The id of the previous halfedge of the triangle for which halfedge with id `eid` is a part.
  if eid %% 3 == 0:
    return eid + 2
  else:
    return eid - 1


iterator iterTriangleEdges*[T](d: Delaunator[T]): tuple[tid: uint32; eid: int32; pid, qid: uint32, p, q: array[2, T]] =
  ## Provides an iterator yielding values for each edge of the triangulation.
  ## The values yielded are the id of the triangle, id of the halfedge chosen
  ## for the edge, id of starting and ending points, an array describing the
  ## point the edge starts at, and an array describing the point the edge ends at.
  var e = 0'i32
  while e < d.triangles.len:
    if e > d.halfedges[e]:
      let
        tid = triangleIdOfEdge(int32(e))
        pid = d.triangles[e]
        qid = d.triangles[nextHalfedge(int32(e))]
        p = [d.coords[2 * pid], d.coords[2 * pid + 1]]
        q = [d.coords[2 * qid], d.coords[2 * qid + 1]]
      yield (tid, e, pid, qid, p, q)
    inc e


func pointIdsOfTriangle*(d: Delaunator, tid: uint32): seq[uint32] =
  ## The point ids composing the triangle with id `tid`.
  return map(halfedgeIdsOfTriangle(tid), proc(eid: int32): uint32 = d.triangles[eid])


iterator iterTriangles*[T](d: Delaunator[T]): tuple[tid: uint32; pid, qid, rid: uint32; p, q, r: array[2, T]] =
  ## Provides an iterator yielding values for each triangle of the triangulation.
  ## The values yielded are the id of the triangle, the ids of the points
  ## comprising the triangle, and three arrays, each describing a point of the triangle.
  var t = 0'u32
  while t < floorDiv(uint32(d.triangles.len), 3'u32):
    let
      pids = pointIdsOfTriangle(d, t)
      p1 = [d.coords[2 * pids[0]], d.coords[2 * pids[0] + 1]]
      p2 = [d.coords[2 * pids[1]], d.coords[2 * pids[1] + 1]]
      p3 = [d.coords[2 * pids[2]], d.coords[2 * pids[2] + 1]]
    yield (t, pids[0], pids[1], pids[2], p1, p2, p3)
    inc t


func triangleIdsAdjacentToTriangle*(d: Delaunator, tid: uint32): seq[uint32] =
  ## The triangle ids adjacent to triangle with id `tid`. A seq of length 2 or 3.
  var tids = collect(newSeqOfCap(3)):
    for eid in halfedgeIdsOfTriangle(tid):
      let opposite = d.halfedges[eid]
      if opposite >= 0: triangleIdOfEdge(opposite)
  return tids


func triangleCentroid*[T](d: Delaunator[T], tid: uint32): array[2, T] =
  ## The centroid of triangle with id `tid`.
  let
    pids = pointIdsOfTriangle(d, tid)
    p1 = [d.coords[2 * pids[0]], d.coords[2 * pids[0] + 1]]
    p2 = [d.coords[2 * pids[1]], d.coords[2 * pids[1] + 1]]
    p3 = [d.coords[2 * pids[2]], d.coords[2 * pids[2] + 1]]
    x = (p1[0] + p2[0] + p3[0]) / 3.0
    y = (p1[1] + p2[1] + p3[1]) / 3.0
  return [x, y]


func polygonCentroid*[T](polygon: seq[array[2, T]]): array[2, T] =
  ## The centroid of points defining `polygon`.
  var
    x, y, k: T = 0.0
    a, b: array[2, T]

  b = polygon[^1]
  for i in 0 ..< polygon.len:
    a = b
    b = polygon[i]
    let c = a[0] * b[1] - b[0] * a[1]
    k += c
    x += (a[0] + b[0]) * c
    y += (a[1] + b[1]) * c

  k *= 3.0
  return [x / k, y / k]


func triangleCircumcenter*[T](d: Delaunator[T], tid: uint32): array[2, T] =
  ## The circumcenter of triangle with id `tid`.
  let
    pids = pointIdsOfTriangle(d, tid)
    p1 = [d.coords[2 * pids[0]], d.coords[2 * pids[0] + 1]]
    p2 = [d.coords[2 * pids[1]], d.coords[2 * pids[1] + 1]]
    p3 = [d.coords[2 * pids[2]], d.coords[2 * pids[2] + 1]]
    (cx, cy) = circumcenter(p1[0], p1[1], p2[0], p2[1], p3[0], p3[1])
  return [cx, cy]


func edgeIdsAroundPoint*(d: Delaunator, eid: int32): seq[int32] =
  ## The ids of all halfedges pointing to the point that halfedge `eid` points to.
  var
    start = eid
    incomming = start
    edgeIds = newSeqOfCap[int32](10) # TODO: Not sure what's a good starting capacity
  while true:
    edgeIds.add(incomming)
    incomming = d.halfedges[nextHalfedge(incomming)]
    if incomming == -1 or incomming == start: break
  return edgeIds


iterator iterVoronoiEdges*[T](d: Delaunator[T]): tuple[eid: int32; p, q: array[2, T]] =
  ## Provides an iterator yielding values for each bisecting voronoi edge of the
  ## triangulation. For finite edges, the values yielded are the id of the
  ## halfedge chosen for the bisected edge, an array describing the circumcenter
  ## of the triangle for which that halfedge is a part, and an array describing
  ## the circumcenter of the adjacent triangle for which that halfedge's
  ## compliment is a part. For infinite edges, the id is of a halfedge on the
  ## hull. The first array describes the circumcenter of the triangle for which
  ## that halfedge is a part, and the second, the point projected by the halfedge
  ## origin's rightmost ray onto the delaunator object's defined bounds.
  # First yield edges of finite regions
  var e = 0'i32
  while e < d.triangles.len:
    if e < d.halfedges[e]: # yield only half the edge, excludes halfedges on hull
      let
        p = triangleCircumcenter[T](d, triangleIdOfEdge(int32(e)))
        q = triangleCircumcenter[T](d, triangleIdOfEdge(d.halfedges[e]))
      yield (e, p, q)
    inc e
  var
    # collect the ids of halfedges on hull
    eids = collect(newSeqOfCap(d.hull.len)):
      for id, val in d.halfedges.pairs:
        if val == -1: int32(id)
  # Second, yield an edge representing rightmost projectedray (from
  # persp. of hull point) of each infinite region on the hull.
  for i in countdown(eids.len - 1, 0):
    let
      e = eids[i]
      p = d.triangles[e]
      incomming = siteToLeftmostHalfedge(d, p)
      edgeIds = edgeIdsAroundPoint(d, incomming)
      triangleIds = map(edgeIds, proc(h: int32): uint32 = triangleIdOfEdge(h))
      vertices = map(triangleIds, proc(t: uint32): array[2, T] = triangleCircumcenter(d, t))
      v = p * 4
      pjctd = project(vertices[^1], [d.vectors[v + 2], d.vectors[v + 3]], d.bounds.minX, d.bounds.minY, d.bounds.maxX, d.bounds.maxY)
    if pjctd.isSome:
      yield (e, vertices[^1], pjctd.get)


# TODO: // degenerate case? (1 valid point: return the box)
#    if (i === 0 && this.delaunay.hull.length === 1) {
#      return [this.xmax, this.ymin, this.xmax, this.ymax, this.xmin, this.ymax, this.xmin, this.ymin];
#    }
iterator iterVoronoiRegions*[T](d: Delaunator[T]): tuple[pid: uint32, verts: seq[array[2, T]]] =
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
        incomming = siteToLeftmostHalfedge(d, p)
        edgeIds = edgeIdsAroundPoint(d, incomming)
        triangleIds = map(edgeIds, proc(h: int32): uint32 = triangleIdOfEdge(h))
        vertices = map(triangleIds, proc(t: uint32): array[2, T] = triangleCircumcenter(d, t))
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


# TODO: // degenerate case (1 valid point: return the box)
#    if (i === 0 && this.delaunay.hull.length === 1) {
#      return [this.xmax, this.ymin, this.xmax, this.ymax, this.xmin, this.ymax, this.xmin, this.ymin];
#    }
proc voronoiRegion*[T](d: Delaunator[T], pid: uint32): tuple[pid: uint32, verts: seq[array[2, T]]] =
  ## Returns the voronoi region of site with id `pid`. The tuple returned
  ## includes the id of the point to which the region belongs, and
  ## a seq of vertices describing the region's polygon. Infinite regions are
  ## clipped to the bounds as set on the delaunator object, defaulting to the
  ## minimum and maximum extents of the points provided. Clipped regions which
  ## are completely out of bounds will yield empty vertices, wheras finite
  ## regions are always yielded whether in bounds or not, and are not clipped.
  let
    incomming = siteToLeftmostHalfedge(d, pid)
    edgeIds = edgeIdsAroundPoint(d, incomming)
    triangleIds = map(edgeIds, proc(h: int32): uint32 = triangleIdOfEdge(h))
    vertices = map(triangleIds, proc(t: uint32): array[2, T] = triangleCircumcenter(d, t))
    inHullAt = find(d.hull, pid)
  if inHullAt == -1: # Not a hull site
    return (pid, vertices)
  else:
    # must clip infinite region
    let
      v = pid * 4 # index into ray vectors
      ply = InfConvexPoly[T](points: vertices, v0: [d.vectors[v], d.vectors[v + 1]], vn: [d.vectors[v + 2], d.vectors[v + 3]])
    return (pid, clipInfinite[T](ply, d.bounds.minX, d.bounds.minY, d.bounds.maxX, d.bounds.maxY))


iterator iterPoints*[T](d: Delaunator[T]): tuple[pid: uint32, p: array[2, T]] =
  ## Provides an iterator yielding values for each point of the triangulation.
  ## The values yielded are the id of the point, and an array describing the
  ## point's location.
  var p = 0'u32
  while p < floorDiv(uint32(d.coords.len), 2'u32):
    yield (p, [d.coords[2 * p], d.coords[2 * p + 1]])
    inc p


iterator iterHullPoints*[T](d: Delaunator[T]): tuple[hid: uint32, pid: uint32, p: array[2, T]] =
  ## Provides an iterator yielding values for each point of the triangulation's
  ## hull. The values yielded are the hull id (index of the point in `d.hull`),
  ## the id of the point, and an array describing the point's location.
  var hid = 0'u32
  while hid < uint32(d.hull.len):
    let pid = d.hull[hid]
    yield (hid, pid, [d.coords[2 * pid], d.coords[2 * pid + 1]])
    inc hid


iterator iterHullEdges*[T](d: Delaunator[T]): tuple[hid: uint32; eid: int32; pid, qid: uint32; p, q: array[2, T]] =
  ## Provides an iterator yielding values for each edge of the triangulation's
  ## hull. The values yielded are the hull id (index of the point in `d.hull` the
  ## edge starts at), the halfedge id, the point id from which the edge starts,
  ## the point id at which the edge ends, an array describing the point the edge
  ## starts at, and an array describing the point the edge ends at.
  var
    # collect the ids of halfedges on hull
    eids = collect(newSeqOfCap(d.hull.len)):
      for id, val in d.halfedges.pairs:
        if val == -1: int32(id)
    hid = 0'u32
    e = eids.len - 1
  while hid < uint32(d.hull.len):
    let
      eid = eids[e]
      pid = d.hull[hid]
      qid = hullNext(d, pid)
      p = [d.coords[2 * pid], d.coords[2 * pid + 1]]
      q = [d.coords[2 * qid], d.coords[2 * qid + 1]]
    yield (hid, eid, pid, qid, p, q)
    inc hid
    dec e
