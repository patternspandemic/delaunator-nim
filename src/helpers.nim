#[
function edgesOfTriangle(t) { return [3 * t, 3 * t + 1, 3 * t + 2]; }
function triangleOfEdge(e)  { return Math.floor(e / 3); }

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
halfedgeIdsOfTriange
triangleIdOfEdge
- nextHalfedge
- prevHalfedge
pointIdsOfTriangle
triangleIdsAdjacentToTriangle
iterPoints
iterHullPoints
iterHullEdges
iterTriangleEdges
iterTriangles
triangleCircumenter
iterVoronoiEdges
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

func nextHalfedge*(e: int32): int32 =
  if e %% 3 == 2:
    result = e - 2
  else:
    result = e + 1


func prevHalfedge*(e: int32): int32 =
  if e %% 3 == 0:
    result = e + 2
  else:
    result = e - 1
