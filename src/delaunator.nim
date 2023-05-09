
import std/math
from std/fenv import epsilon


var EDGE_STACK: array[512, uint32]

type
  Delaunator*[T] = ref object
    coords*: seq[T]
    triangles*: seq[uint32] # trimmed version of d_triangles
    halfedges*: seq[int32]  # trimmed version of d_halfedges
    hull*: seq[uint32]

    # Arrays that will store the triangulation graph
    trianglesLen: int32
    d_triangles: seq[uint32]
    d_halfedges: seq[int32]

    # Temporary arrays for tracking the edges of the advancing convex hull
    d_hashSize: int
    d_hullStart: int
    d_hullPrev: seq[uint32]  # edge to prev edge
    d_hullNext: seq[uint32]  # edge to next edge
    d_hullTri:  seq[uint32]  # edge to adjacent triangle
    d_hullHash: seq[int32]   # angular edge hash

    # Temporary arrays for sorting points
    d_ids:   seq[uint32]
    d_dists: seq[T]


func defaultGetX[P, T](p: P): T = p[0]
func defaultGetY[P, T](p: P): T = p[1]


proc swap(arr: var seq[uint32]; i, j: int) =
  let tmp = arr[i]
  arr[i] = arr[j]
  arr[j] = tmp


func pseudoAngle[F](dx, dy: F): F =
  let p = dx / (dx.abs + dy.abs)
  result = (if dy > 0: 3 - p else: 1 + p)


func dist[F](ax, ay, bx, by: F): F =
  let
    dx = ax - bx
    dy = ay - by
  result = dx * dx + dy * dy


func inCircle[F](ax, ay, bx, by, cx, cy, px, py: F): bool =
  let
    dx = ax - px
    dy = ay - py
    ex = bx - px
    ey = by - py
    fx = cx - px
    fy = cy - py

    ap = dx * dx + dy * dy
    bp = ex * ex + ey * ey
    cp = fx * fx + fy * fy

    result = dx * (ey * cp - bp * fy) -
             dy * (ex * cp - bp * fx) +
             ap * (ex * fy - ey * fx) < 0


func circumradius[F](ax, ay, bx, by, cx, cy: F): F =
  let
    dx = bx - ax
    dy = by - ay
    ex = cx - ax
    ey = cy - ay

    bl = dx * dx + dy * dy
    cl = ex * ex + ey * ey
    d = 0.5 / (dx * ey - dy * ex)

    x = (ey * bl - dy * cl) * d
    y = (dx * cl - ex * bl) * d

  result = x * x + y * y


func circumcenter[F](ax, ay, bx, by, cx, cy: F): tuple[x, y: F] =
  let
    dx = bx - ax
    dy = by - ay
    ex = cx - ax
    ey = cy - ay

    bl = dx * dx + dy * dy
    cl = ex * ex + ey * ey
    d = 0.5 / (dx * ey - dy * ex)

    x = ax + (ey * bl - dy * cl) * d
    y = ay + (dx * cl - ex * bl) * d

  result = (x, y)


proc quicksort[F](ids: var seq[uint32]; dists: seq[F]; left, right: int) =
  if (right - left <= 20):
    var i = left + 1
    while i <= right:
      let
        temp = ids[i]
        tempDist = dists[temp]
      var j = i - 1
      while j >= left and dists[ids[j]] > tempDist:
        ids[j + 1] = ids[j]
        dec j
      ids[j + 1] = temp
      inc i
  else:
    let median = ashr(left + right, 1)
    var
      i = left + 1
      j = right
    swap(ids, median, i)
    if dists[ids[left]] > dists[ids[right]]: swap(ids, left, right)
    if dists[ids[i]] > dists[ids[right]]: swap(ids, i, right)
    if dists[ids[left]] > dists[ids[i]]: swap(ids, left, i)

    let
      temp = ids[i]
      tempDist = dists[temp]
    while true:
      while true:
        inc i
        if dists[ids[i]] < tempDist: break
      while true:
        dec j
        if dists[ids[j]] > tempDist: break
      if j < i: break
      swap(ids, i, j)
    ids[left + 1] = ids[j]
    ids[j] = temp

    if right - i + 1 >= j - left:
      quicksort(ids, dists, i, right)
      quicksort(ids, dists, left, j - 1)
    else:
      quicksort(ids, dists, left, j - 1)
      quicksort(ids, dists, i, right)


proc update(this: var Delaunator) =

  # TODO, move this proc down to where cx, cy are delcared.
  let
    (u_cx, u_cy) = (1.0, 2.0) # circumcenter(....)
  proc u_hashKey(x, y: SomeFloat): int32 =
    return floor(pseudoAngle(x - u_cx, y - u_cy) * this.d_hashSize) %% this.d_hashSize


  proc u_link(a, b: int32) =
    this.d_halfedges[a] = b
    if b != -1: this.d_halfedges[b] = a


  # Lots of int type casting due to disconnect between indexing and unsigned storage types used for halfedges/hull checks.
  proc u_legalize(a: var int32): uint32 =
    var
      triangles = this.d_triangles
      halfedges = this.d_halfedges
      coords = this.coords
      i = 0
      ar = 0'i32

    while true:
      block outerWhile:
        let b = halfedges[a]

        #[
        if the pair of triangles doesn't satisfy the Delaunay condition
        (p1 is inside the circumcircle of [p0, pl, pr]), flip them,
        then do the same check/flip recursively for the new pair of triangles

                pl                    pl
               /||\                  /  \
            al/ || \bl            al/    \a
             /  ||  \              /      \
            /  a||b  \    flip    /___ar___\
          p0\   ||   /p1   =>   p0\---bl---/p1
             \  ||  /              \      /
            ar\ || /br             b\    /br
               \||/                  \  /
                pr                    pr

        ]#

        let a0 = a - a %% 3
        ar = a0 + (a + 2) %% 3

        if b == -1: # convex hull edge
          if i == 0: break
          dec i
          a = int32(EDGE_STACK[i])
          continue

        let
          b0 = b - b %% 3
          al = a0 + (a + 1) %% 3
          bl = b0 + (b + 2) %% 3

          p0 = triangles[ar]
          pr = triangles[a]
          pl = triangles[al]
          p1 = triangles[bl]

          illegal = inCircle(
            coords[2 * p0], coords[2 * p0 + 1],
            coords[2 * pr], coords[2 * pr + 1],
            coords[2 * pl], coords[2 * pl + 1],
            coords[2 * p1], coords[2 * p1 + 1])

        if illegal:
          triangles[a] = p1
          triangles[b] = p0

          let hbl = halfedges[bl]

          # edge swapped on the other side of the hull (rare); fix the halfedge reference
          if hbl == -1:
            var e = this.d_hullStart
            while true:
              if int32(this.d_hullTri[e]) == bl:
                this.d_hullTri[e] = uint32(a)
                break outerWhile
              e = int32(this.d_hullPrev[e])
              if e == this.d_hullStart: break
          u_link(a, hbl)
          u_link(b, halfedges[ar])
          u_link(ar, bl)

          let br = b0 + (b + 1) %% 3

          # don't worry about hitting the cap: it can only happen on extremely degenerate input
          if i < EDGE_STACK.len:
            EDGE_STACK[i] = uint32(br)
            inc i

        else:
          if i == 0: break
          dec i
          a = int32(EDGE_STACK[i])

    return uint32(ar)


  proc u_addTriangle(i0, i1, i2, a, b, c: int32): int32 =
    let t = this.trianglesLen

    this.d_triangles[t] = uint32(i0)
    this.d_triangles[t + 1] = uint32(i1)
    this.d_triangles[t + 2] = uint32(i2)

    u_link(t, a)
    u_link(t + 1, b)
    u_link(t + 2, c)

    this.trianglesLen += 3

    return t


  # TODO: main update code here
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
  result.d_hullStart = 0
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
