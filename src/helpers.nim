# TODO: inline things
#       iterators sans ids?

import std/sequtils, std/sugar

#[

function forEachTriangleEdge(points, delaunay, callback) {
    for (let e = 0; e < delaunay.triangles.length; e++) {
        if (e > delaunay.halfedges[e]) {
            const p = points[delaunay.triangles[e]];
            const q = points[delaunay.triangles[nextHalfedge(e)]];
            callback(e, p, q);
        }
    }
}

function pointsOfTriangle(delaunay, t) {
    return edgesOfTriangle(t)
        .map(e => delaunay.triangles[e]);
}

function forEachTriangle(points, delaunay, callback) {
    for (let t = 0; t < delaunay.triangles.length / 3; t++) {
        callback(t, pointsOfTriangle(delaunay, t).map(p => points[p]));
    }
}

function trianglesAdjacentToTriangle(delaunay, t) {
    const adjacentTriangles = [];
    for (const e of edgesOfTriangle(t)) {
        const opposite = delaunay.halfedges[e];
        if (opposite >= 0) {
            adjacentTriangles.push(triangleOfEdge(opposite));
        }
    }
    return adjacentTriangles;
}

function circumcenter(a, b, c) {
    const ad = a[0] * a[0] + a[1] * a[1];
    const bd = b[0] * b[0] + b[1] * b[1];
    const cd = c[0] * c[0] + c[1] * c[1];
    const D = 2 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]));
    return [
        1 / D * (ad * (b[1] - c[1]) + bd * (c[1] - a[1]) + cd * (a[1] - b[1])),
        1 / D * (ad * (c[0] - b[0]) + bd * (a[0] - c[0]) + cd * (b[0] - a[0])),
    ];
}

function triangleCenter(points, delaunay, t) {
    const vertices = pointsOfTriangle(delaunay, t).map(p => points[p]);
    return circumcenter(vertices[0], vertices[1], vertices[2]);
}

function forEachVoronoiEdge(points, delaunay, callback) {
    for (let e = 0; e < delaunay.triangles.length; e++) {
        if (e < delaunay.halfedges[e]) {
            const p = triangleCenter(points, delaunay, triangleOfEdge(e));
            const q = triangleCenter(points, delaunay, triangleOfEdge(delaunay.halfedges[e]));
            callback(e, p, q);
        }
    }
}

function edgesAroundPoint(delaunay, start) {
    const result = [];
    let incoming = start;
    do {
        result.push(incoming);
        const outgoing = nextHalfedge(incoming);
        incoming = delaunay.halfedges[outgoing];
    } while (incoming !== -1 && incoming !== start);
    return result;
}

function forEachVoronoiCell(points, delaunay, callback) {
    const index = new Map(); // point id to half-edge id
    for (let e = 0; e < delaunay.triangles.length; e++) {
        const endpoint = delaunay.triangles[nextHalfedge(e)];
        if (!index.has(endpoint) || delaunay.halfedges[e] === -1) {
            index.set(endpoint, e);
        }
    }
    for (let p = 0; p < points.length; p++) {
        const incoming = index.get(p);
        const edges = edgesAroundPoint(delaunay, incoming);
        const triangles = edges.map(triangleOfEdge);
        const vertices = triangles.map(t => triangleCenter(points, delaunay, t));
        callback(p, vertices);
    }
}
]#

#[
bounds
- halfedgeIdsOfTriange
- triangleIdOfEdge
- nextHalfedge
- prevHalfedge
- pointIdsOfTriangle
- triangleIdsAdjacentToTriangle
iterPoints
iterHullPoints
iterHullEdges
- iterTriangleEdges
- iterTriangles
- triangleCircumenter
- iterVoronoiEdges
edgeIdsAroundPoint
iterVoronoiRegions
voronoiRegion

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


func halfedgeIdsOfTriangle(t: int32): array[3, int32] =
  ## The halfedge ids of a triangle with id `t`.
  return [3 * t, 3 * t + 1, 3 * t + 2]


func triangleIdOfEdge(e: int32): int32 =
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


iterator iterTriangleEdges*(d: Delaunator): tuple[e: int32, p: array[2, SomeFloat], q: array[2, SomeFloat]] =
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
        p: array[2, SomeFloat] = d.coords[(2 * pid)..(2 * pid + 2)]
        q: array[2, SomeFloat] = d.coords[(2 * qid)..(2 * qid + 2)]
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
      p1: array[2, SomeFloat] = d.coords[(2 * pids[0])..(2 * pids[0] + 2)]
      p2: array[2, SomeFloat] = d.coords[(2 * pids[1])..(2 * pids[1] + 2)]
      p3: array[2, SomeFloat] = d.coords[(2 * pids[2])..(2 * pids[2] + 2)]
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
    p1: array[2, SomeFloat] = d.coords[(2 * pids[0])..(2 * pids[0] + 2)]
    p2: array[2, SomeFloat] = d.coords[(2 * pids[1])..(2 * pids[1] + 2)]
    p3: array[2, SomeFloat] = d.coords[(2 * pids[2])..(2 * pids[2] + 2)]
    (cx, cy) = circumcenter(p1[0], p1[1], p2[0], p2[1], p3[0], p3[1])
  return [cx, cy]

iterator iterVoronoiEdges*(d: Delaunator): tuple[e: int32, p, q: array[2, SomeFloat]] =
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
