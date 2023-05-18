
# Up to date with mapbox/Delaunator at 103acb4564a36ad2dff11dc0135a348f4e8fc149 May 27, 2032

#TODO: inline things

import std/[math, tables]
from std/fenv import epsilon
from std/algorithm import fill

from orient2d import orient2d


var EDGE_STACK: array[512, uint32]

type
  Delaunator*[T] = ref object
    coords*: seq[T]
    minX*, minY*, maxX*, maxY*: T # bounds of coords
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

    # For fast lookup of point id to leftmost imcoming halfedge id
    # Useful for retrieval of adhoc voronoi regions.
    d_pointToLeftmostHalfedgeIndex: Table[uint32, int32]


proc swap(arr: var seq[uint32]; i, j: int) =
  let tmp = arr[i]
  arr[i] = arr[j]
  arr[j] = tmp


# monotonically increases with real angle, but doesn't need expensive trigonometry
func pseudoAngle[F](dx, dy: F): F =
  let p = dx / (dx.abs + dy.abs)
  if dy > 0.0:
    result = (3.0 - p) / 4.0
  else:
    result = (1.0 + p) / 4.0


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

    x: F = ax + (ey * bl - dy * cl) * d
    y: F = ay + (dx * cl - ex * bl) * d

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
        if not (dists[ids[i]] < tempDist): break
      while true:
        dec j
        if not (dists[ids[j]] > tempDist): break
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


include helpers # some depend upon above


proc update[T](this: var Delaunator) =
  # Inner procs passed var 'this' param as 'uthis' due to the need to mutate it.
  # Simply closing over it results in 'cannot be captured / memory safety error'.

  proc u_link(uthis: var Delaunator; a, b: int32) =
    uthis.d_halfedges[a] = b
    if b != -1: uthis.d_halfedges[b] = a


  # Lots of int type casting due to disconnect between indexing and unsigned storage types used for halfedges/hull checks.
  proc u_legalize(uthis: var Delaunator; a: var int32): uint32 =
    var
      i = 0
      ar = 0'i32

    block outerWhile: # TODO: remove block
      while true:
        let b = uthis.d_halfedges[a]

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

        let a0 = a - (a mod 3)
        ar = a0 + ((a + 2) mod 3)

        if b == -1: # convex hull edge
          if i == 0: break
          dec i
          a = int32(EDGE_STACK[i])
          continue

        let
          b0 = b - (b mod 3)
          al = a0 + ((a + 1) mod 3)
          bl = b0 + ((b + 2) mod 3)

          p0 = uthis.d_triangles[ar]
          pr = uthis.d_triangles[a]
          pl = uthis.d_triangles[al]
          p1 = uthis.d_triangles[bl]

          illegal = inCircle(
            uthis.coords[2 * p0], uthis.coords[2 * p0 + 1],
            uthis.coords[2 * pr], uthis.coords[2 * pr + 1],
            uthis.coords[2 * pl], uthis.coords[2 * pl + 1],
            uthis.coords[2 * p1], uthis.coords[2 * p1 + 1])

        if illegal:
          uthis.d_triangles[a] = p1
          uthis.d_triangles[b] = p0

          let hbl = uthis.d_halfedges[bl]

          # edge swapped on the other side of the hull (rare); fix the halfedge reference
          if hbl == -1:
            var e = uthis.d_hullStart
            while true: #FIXME: issue13, issue44 hang here
              if int32(uthis.d_hullTri[e]) == bl:
                uthis.d_hullTri[e] = uint32(a)
                break #outerWhile
              e = int32(uthis.d_hullPrev[e])
              if e == uthis.d_hullStart: break
          u_link(uthis, a, hbl)
          u_link(uthis, b, uthis.d_halfedges[ar])
          u_link(uthis, ar, bl)

          let br = b0 + ((b + 1) mod 3)

          # don't worry about hitting the cap: it can only happen on extremely degenerate input
          if i < EDGE_STACK.len:
            EDGE_STACK[i] = uint32(br)
            inc i

        else:
          if i == 0: break
          dec i
          a = int32(EDGE_STACK[i])

    return uint32(ar)


  proc u_addTriangle(uthis: var Delaunator; i0, i1, i2, a, b, c: int): int32 =
    let t = uthis.trianglesLen

    uthis.d_triangles[t] = uint32(i0)
    uthis.d_triangles[t + 1] = uint32(i1)
    uthis.d_triangles[t + 2] = uint32(i2)

    u_link(uthis, t, int32(a))
    u_link(uthis, t + 1, int32(b))
    u_link(uthis, t + 2, int32(c))

    uthis.trianglesLen += 3

    return t


  # The main update code begins here (nested procs mostly above)
  let n = ashr(this.coords.len, 1)

  # populate an array of point indices; calculate input data bbox
  this.minX = Inf
  this.minY = Inf
  this.maxX = NegInf
  this.maxY = NegInf

  for i in 0 ..< n:
    let
      x = this.coords[2 * i]
      y = this.coords[2 * i + 1]
    if x < this.minX: this.minX = x
    if y < this.minY: this.minY = y
    if x > this.maxX: this.maxX = x
    if y > this.maxY: this.maxY = y
    this.d_ids[i] = uint32(i)
  let
    cx = (this.minX + this.maxX) / 2
    cy = (this.minY + this.maxY) / 2

  var
    i0, i1, i2: int
    i0x, i0y, i1x, i1y, i2x, i2y = T(0) # Temp init to something for case of empty coords.

  # pick a seed point close to the center
  var minDist = Inf
  for i in 0 ..< n:
    let d = dist(cx, cy, this.coords[2 * i], this.coords[2 * i + 1])
    if d < minDist:
      i0 = i
      minDist = d

  if this.coords.len > 0:
    i0x = this.coords[2 * i0]
    i0y = this.coords[2 * i0 + 1]

  # find the point closest to the seed
  minDist = Inf
  for i in 0 ..< n:
    if i == i0: continue
    let d = dist(i0x, i0y, this.coords[2 * i], this.coords[2 * i + 1])
    if d < minDist and d > 0:
      i1 = i
      minDist = d

  if this.coords.len > 0:
    i1x = this.coords[2 * i1]
    i1y = this.coords[2 * i1 + 1]

  var minRadius = Inf

  # find the third point which forms the smallest circumcircle with the first two
  for i in 0 ..< n:
    if i == i0 or i == i1: continue
    let r = circumradius(i0x, i0y, i1x, i1y, this.coords[2 * i], this.coords[2 * i + 1])
    if r < minRadius:
      i2 = i
      minRadius = r

  if this.coords.len > 0:
    i2x = this.coords[2 * i2]
    i2y = this.coords[2 * i2 + 1]

  if minRadius == Inf:
    # order collinear points by dx (or dy if all x are identical)
    # and return the list as a hull
    for i in 0 ..< n:
      let
        xcrd = this.coords[2 * i] - this.coords[0]
        ycrd = this.coords[2 * i + 1] - this.coords[1]
      if xcrd != 0:
        this.d_dists[i] = xcrd
      else:
        this.d_dists[i] = ycrd
    quicksort(this.d_ids, this.d_dists, 0, n - 1)
    var hull = newSeq[uint32](n)
    var
      j = 0
      d0 = NegInf
    for i in 0 ..< n:
      let
        id = this.d_ids[i]
        d = this.d_dists[id]
      if d > d0:
        hull[j] = uint32(id)
        inc j
        d0 = d
    this.hull = hull[0 ..< j]
    this.triangles = newSeqOfCap[uint32](0)
    this.halfedges = newSeqOfCap[int32](0)
    return

  # swap the order of the seed points for counter-clockwise orientation
  if orient2d(i0x, i0y, i1x, i1y, i2x, i2y) < 0:
    let
      i = i1
      x = i1x
      y = i1y
    i1 = i2
    i1x = i2x
    i1y = i2y
    i2 = i
    i2x = x
    i2y = y

  let
    (u_cx, u_cy) = circumcenter[T](i0x, i0y, i1x, i1y, i2x, i2y)

  # defined here to close over u_cx, u_cy
  proc u_hashKey(uthis: var Delaunator; x, y: SomeFloat): int32 =
    return int32(floor(pseudoAngle(x - u_cx, y - u_cy) * float(uthis.d_hashSize)) mod float(uthis.d_hashSize))

  for i in 0 ..< n:
    this.d_dists[i] = dist(this.coords[2 * i], this.coords[2 * i + 1], u_cx, u_cy)

  # sort the points by distance from the seed triangle circumcenter
  quicksort(this.d_ids, this.d_dists, 0, n - 1)

  # set up the seed triangle as the starting hull
  this.d_hullStart = i0
  var hullSize = 3

  this.d_hullNext[i0] = uint32(i1); this.d_hullPrev[i2] = uint32(i1)
  this.d_hullNext[i1] = uint32(i2); this.d_hullPrev[i0] = uint32(i2)
  this.d_hullNext[i2] = uint32(i0); this.d_hullPrev[i1] = uint32(i0)

  this.d_hullTri[i0] = 0
  this.d_hullTri[i1] = 1
  this.d_hullTri[i2] = 2

  this.d_hullHash.fill(-1)
  this.d_hullHash[u_hashKey(this, i0x, i0y)] = int32(i0)
  this.d_hullHash[u_hashKey(this, i1x, i1y)] = int32(i1)
  this.d_hullHash[u_hashKey(this, i2x, i2y)] = int32(i2)

  #this.trianglesLen = 0 # already init'd in fromCoords
  discard u_addTriangle(this, i0, i1, i2, -1, -1, -1)

  var
    xp, yp = this.coords[0] # just so xp & yp are of correct type
  for k in 0 ..< this.d_ids.len:
    let
      i = this.d_ids[k]
      x = this.coords[2 * i]
      y = this.coords[2 * i + 1]

    # skip near-duplicate points
    if k > 0 and abs(x - xp) <= epsilon(float64) and abs(y - yp) <= epsilon(float64): continue
    xp = x
    yp = y

    # skip seed triangle points
    if i == uint32(i0) or i == uint32(i1) or i == uint(i2): continue

    # find a visible edge on the convex hull using edge hash
    var
      start = 0
      key = u_hashKey(this, x, y)
    for j in 0 ..< this.d_hashSize:
      start = this.d_hullHash[(key + j) mod this.d_hashSize]
      if start != -1 and uint32(start) != this.d_hullNext[start]: break

    start = int(this.d_hullPrev[start])
    var
      e = start
      q = int(this.d_hullNext[e])
    while orient2d(x, y, this.coords[2 * e], this.coords[2 * e + 1], this.coords[2 * q], this.coords[2 * q + 1]) >= 0:
      e = q
      if e == start:
        e = -1
        break
      q = int(this.d_hullNext[e])
    if e == -1: continue # likely a near-duplicate point; skip it

    # add the first triangle from the point
    var t = u_addTriangle(this, e, int(i), int(this.d_hullNext[e]), -1, -1, int(this.d_hullTri[e]))

    # recursively flip triangles from the point until they satisfy the Delaunay condition
    var tmp = t + 2 # use mutable tmp to make u_legalize happy
    this.d_hullTri[i] = u_legalize(this, tmp)
    this.d_hullTri[e] = uint32(t) # keep track of boundary triangles on the hull
    inc hullSize

    # walk forward through the hull, adding more triangles and flipping recursively
    var n = int(this.d_hullNext[e])
    q = int(this.d_hullNext[n])
    while orient2d(x, y, this.coords[2 * n], this.coords[2 * n + 1], this.coords[2 * q], this.coords[2 * q + 1]) < 0:
      t = u_addTriangle(this, n, int(i), q, int(this.d_hullTri[i]), -1, int(this.d_hullTri[n]))
      tmp = t + 2 # use mutable tmp to make u_legalize happy
      this.d_hullTri[i] = u_legalize(this, tmp)
      this.d_hullNext[n] = uint32(n) # mark as removed
      dec hullSize
      n = q
      q = int(this.d_hullNext[n])

    # walk backward from the other side, adding more triangles and flipping
    if e == start:
      q = int(this.d_hullPrev[e])
      while orient2d(x, y, this.coords[2 * q], this.coords[2 * q + 1], this.coords[2 * e], this.coords[2 * e + 1]) < 0:
        t = u_addTriangle(this, q, int(i), e, -1, int(this.d_hullTri[e]), int(this.d_hullTri[q]))
        tmp = t + 2 # use mutable tmp to make u_legalize happy
        discard u_legalize(this, tmp)
        this.d_hullTri[q] = uint32(t)
        this.d_hullNext[e] = uint32(e) # mark as removed
        dec hullSize
        e = q
        q = int(this.d_hullPrev[e])

    # update the hull indices
    this.d_hullPrev[i] = uint32(e)
    this.d_hullStart = e
    this.d_hullPrev[n] = uint32(i)
    this.d_hullNext[e] = uint32(i)
    this.d_hullNext[i] = uint32(n)

    # save the two new edges in the hash table
    this.d_hullHash[u_hashKey(this, x, y)] = int32(i)
    this.d_hullHash[u_hashKey(this, this.coords[2 * e], this.coords[2 * e + 1])] = int32(e)

  this.hull = newSeq[uint32](hullSize)
  var e = this.d_hullStart
  for i in 0 ..< hullSize:
    this.hull[i] = uint32(e)
    e = int(this.d_hullNext[e])

  # Build the index of point id to leftmost incoming halfedge
  clear(this.d_pointToLeftmostHalfedgeIndex)
  var he: int32 = 0
  while true:
    let endpoint = this.d_triangles[nextHalfedge(he)]
    if not hasKey(this.d_pointToLeftmostHalfedgeIndex, endpoint) or this.d_halfedges[he] == -1:
      this.d_pointToLeftmostHalfedgeIndex[endpoint] = he
    inc he
    if not he < this.trianglesLen: break

  # trim typed triangle mesh arrays
  this.triangles = this.d_triangles[0 ..< this.trianglesLen]
  this.halfedges = this.d_halfedges[0 ..< this.trianglesLen]


proc fromCoords*[T](coordinates: seq[T]): Delaunator[T] =
  # Could not figure out how to constrain T to SomeFloat, so using cound to test coordinates as
  # a seq[float32|float64]. Errors for example when coordinates is seq[int32].
  when not compiles(count[T](coordinates, 0.0)):
    {.error: "Coordinates must be seq[float32] or seq[float64] but got " & $typeof(coordinates) .}
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
  update[T](result)


func defaultGetX[P, T](p: P): T =
  ## Default getX proc for `fromPoints`. Coerces to `T`.
  T(p[0])


func defaultGetY[P, T](p: P): T =
  ## Default getY proc for `fromPoints`. Coerces to `T`.
  T(p[1])


proc fromPoints*[P, T](points: seq[P]; getX: proc (p: P): T = defaultGetX; getY: proc (p: P): T = defaultGetY): Delaunator[T] =
  var
    coords = newSeq[T](points.len * 2)
  for i, point in points:
      coords[2 * i] = getX(point)
      coords[2 * i + 1] = getY(point)
  fromCoords[T](coords)



when isMainModule:
  #[
  type
    Point32 = tuple[x, y: float32]
    Point64 = tuple[x, y: float64]

  var
    d32: Delaunator[float32]
    points32: seq[Point32]
    d64: Delaunator[float64]
    points64: seq[Point64]

  points32 = @[(1.0'f32, 1.0'f32), (-1.0'f32, 1.0'f32), (0.2'f32, 0.1'f32), (0.0'f32, -1.0'f32)]
  points64 = @[(1.0'f64, 1.0'f64), (-1.0'f64, 1.0'f64), (0.2'f64, 0.1'f64), (0.0'f64, -1.0'f64)]
  d32 = fromPoints[Point32, float32](points32)
  d64 = fromPoints[Point64, float64](points64)

  echo repr(d32)
  echo ""
  echo "Coords: :", d32.coords
  echo "Triangles: ", d32.triangles
  echo "Halfedges: ", d32.halfedges
  echo "Hull: ", d32.hull
  echo " - - - - - - - - "
  echo repr(d64)
  echo ""
  echo "Coords: :", d64.coords
  echo "Triangles: ", d64.triangles
  echo "Halfedges: ", d64.halfedges
  echo "Hull: ", d64.hull
  ]#

  #include "../tests/fixtures/ukraine"
  #include "../tests/fixtures/issue44"
  include "../tests/fixtures/robustness4"

  var d = fromPoints[array[2, float64], float64](robustness4)
  #echo "Coords type: ", typeof(d.coords)
  echo repr(d)
  # echo "Coords: :", d.coords
  # echo "Triangles: ", d.triangles
  # echo "Halfedges: ", d.halfedges
  # echo "Hull: ", d.hull
