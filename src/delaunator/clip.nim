
# Clipping code based on https://observablehq.com/@mbostock/to-infinity-and-back-again
# as well as d3-delaunay.

import std/[options, sequtils]


type
  InfConvexPoly*[T] = object
    points*: seq[array[2, T]]
    v0*, vn*: array[2, T]


# TODO: Find where this is called. in linked example
#[
function centroid(polygon) {
  var i = -1,
      n = polygon.length,
      x = 0,
      y = 0,
      a,
      b = polygon[n - 1],
      c,
      k = 0;

  while (++i < n) {
    a = b;
    b = polygon[i];
    k += c = a[0] * b[1] - b[0] * a[1];
    x += (a[0] + b[0]) * c;
    y += (a[1] + b[1]) * c;
  }

  return k *= 3, [x / k, y / k];
}
]#
func centroid[T](polygon: seq[array[2, T]]): array[2, T] =
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


func project*[T](p, v: array[2, T], xMin, yMin, xMax, yMax: T): Option[array[2, T]] =
  var
    t = Inf
    c, x, y: T

  if v[1] < 0.0: # top
    if p[1] <= yMin: return # none(array[2,T])
    c = (yMin - p[1]) / v[1]
    if c < t:
      y = yMin
      x = p[0] + c * v[0]
      t = c
  elif v[1] > 0.0: # bottom
    if p[1] >= yMax: return  # none(array[2,T])
    c = (yMax - p[1]) / v[1]
    if c < t:
      y = yMax
      x = p[0] + c * v[0]
      t = c

  if v[0] > 0.0: # right
    if p[0] >= xMax: return  # none(array[2,T])
    c = (xMax - p[0]) / v[0]
    if c < t:
      x = xMax
      y = p[1] + c * v[1]
      t = c
  elif v[0] < 0.0: # left
    if p[0] <= xMin: return # none(array[2,T])
    c = (xMin - p[0]) / v[0]
    if c < t:
      x = xMin
      y = p[1] + c * v[1]
      t = c

  return some([x, y])


proc clipInfinite*[T](polygon: InfConvexPoly[T], xMin, yMin, xMax, yMax: T): seq[array[2, T]] =
  let # so inner procs can close over bounds
    xMin = xMin
    yMin = yMin
    xMax = xMax
    yMax = yMax


  func clockwise(p, q, r: array[2, T]): bool =
    return (q[0] - p[0]) * (r[1] - p[1]) < (q[1] - p[1]) * (r[0] - p[0])


  func contains(poly: InfConvexPoly[T], p: array[2, T]): bool =
    var
      p0, p1: array[2, T]
    p1 = poly.points[0]
    if clockwise(p, [p1[0] + poly.v0[0], p1[1] + poly.v0[1]], p1): return false
    for i in 1 ..< poly.points.len:
      p0 = p1
      p1 = poly.points[i]
      if clockwise(p, p0, p1): return false
    if clockwise(p, p1, [p1[0] + poly.vn[0], p1[1] + poly.vn[1]]): return false
    return true


  let
    insideTop = proc (p: array[2, T]): bool = p[1] > yMin
    insideRight = proc (p: array[2, T]): bool = p[0] < xMax
    insideBottom = proc (p: array[2, T]): bool = p[1] < yMax
    insideLeft = proc (p: array[2, T]): bool = p[0] > xMin

    intersectTop = proc (p, q: array[2, T]): array[2, T] = [p[0] + ( q[0] - p[0] ) * ( yMin - p[1] ) / ( q[1] - p[1] ), yMin]
    intersectRight = proc (p, q: array[2, T]): array[2, T] = [xMax, p[1] + ( q[1] - p[1] ) * ( xMax - p[0] ) / ( q[0] - p[0] )]
    intersectBottom = proc (p, q: array[2, T]): array[2, T] = [p[0] + ( q[0] - p[0] ) * ( yMax - p[1] ) / ( q[1] - p[1] ), yMax]
    intersectLeft = proc (p, q: array[2, T]): array[2, T] = [xMin, p[1] + ( q[1] - p[1] ) * ( xMin - p[0] ) / ( q[0] - p[0] )]


  proc clipper(
    inside: proc (p: array[2, T]): bool,
    intersect: proc (p, q: array[2, T]): array[2, T]
       ): proc (poly: seq[array[2, T]]): seq[array[2, T]] =
    return proc (subject: seq[array[2, T]]): seq[array[2, T]] =
      var p = newSeqOfCap[array[2, T]](subject.len + 5) # FIXME: Hmm, At most adding 5 points to poly durring clip?
      if subject.len == 0: return p
      var
        p0, p1: array[2, T]
        t0, t1: bool
      p1 = subject[^1]
      t1 = inside(p1)
      for i in 0 ..< subject.len:
        p0 = p1
        p1 = subject[i]
        t0 = t1
        t1 = inside(p1)
        if t1 != t0: p.add(intersect(p0, p1))
        if t1: p.add(p1)
      return p


  let
    top = clipper(insideTop, intersectTop)
    right = clipper(insideRight, intersectRight)
    bottom = clipper(insideBottom, intersectBottom)
    left = clipper(insideLeft, intersectLeft)


  proc clip(subject: seq[array[2, T]]): seq[array[2, T]] =
    return left(bottom(right(top(subject))))


  func sidecode(p: array[2, T]): int8 =
    return (if p[0] == xMin: int8(1) else: (if p[0] == xMax: int8(2) else: int8(0))) or (if p[1] == yMin: int8(4) else: (if p[1] == yMax: int8(8) else: int8(0)))


  var
    ply = polygon.points
    p: Option[array[2, T]]
    n: int
  p = project[T](ply[0], polygon.v0, xMin, yMin, xMax, yMax)
  if p.isSome: ply.insert(@[p.get]) # FIXME: Better way than insert?
  p = project[T](ply[^1], polygon.vn, xMin, yMin, xMax, yMax)
  if p.isSome: ply.insert(@[p.get])
  ply = clip(ply)
  n = ply.len
  if n > 0:
    var
      i = 0
      c0, c1: int8
    c1 = sidecode(ply[n - 1])
    while i < n:
      c0 = c1
      c1 = sidecode(ply[i])
      if (c0 != 0) and (c1 != 0):
        while c0 != c1:
          var c: array[2, T]
          case c0
          of 0b0101: # top-left
            c0 = 0b0100
            continue
          of 0b0100: # top
            c0 = 0b0110
            c = [xMax, yMin]
          of 0b0110: # top-right
            c0 = 0b0010
            continue
          of 0b0010: # right
            c0 = 0b1010
            c = [xMax, yMax]
          of 0b1010: # bottom-right
            c0 = 0b1000
            continue
          of 0b1000: # bottom
            c0 = 0b1001
            c = [xMin, yMax]
          of 0b1001: # bottom-left
            c0 = 0b0001
            continue
          of 0b0001: # left
            c0 = 0b0101
            c = [xMin, yMin]
          else:
            discard
          if contains(polygon, c):
            ply.insert(@[c], i)
            inc n
            inc i
      inc i

  elif contains(polygon, [(xMin + xMax) / 2.0, (yMin + yMax) / 2.0]):
    ply.insert(@[[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax]])

  return ply
