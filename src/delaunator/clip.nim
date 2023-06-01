
# Clipping code based on https://observablehq.com/@mbostock/to-infinity-and-back-again

type
  InfConvexPoly[T] = object
    points: seq[array[2, T]]
    v0, vn: array[2, T]


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


proc clipInfinite*[T](polygon: InfConvexPoly[T], xMin, yMin, xMax, yMax: T): seq[array[2, T]] =
  #[
  function project([x0, y0], [vx, vy]) {
    let t = Infinity, c, x, y;
    if (vy < 0) { // top
      if (y0 <= ymin) return;
      if ((c = (ymin - y0) / vy) < t) y = ymin, x = x0 + (t = c) * vx;
    } else if (vy > 0) { // bottom
      if (y0 >= ymax) return;
      if ((c = (ymax - y0) / vy) < t) y = ymax, x = x0 + (t = c) * vx;
    }
    if (vx > 0) { // right
      if (x0 >= xmax) return;
      if ((c = (xmax - x0) / vx) < t) x = xmax, y = y0 + (t = c) * vy;
    } else if (vx < 0) { // left
      if (x0 <= xmin) return;
      if ((c = (xmin - x0) / vx) < t) x = xmin, y = y0 + (t = c) * vy;
    }
    return [x, y];
  }
  ]#
  func project(p, v: array[2, T]): array[2, T] =
    var
      t = Inf
      c, x, y: T

    if v[1] < 0.0: # top
      if p[1] <= yMin: return nil
      c = (yMin - p[1]) / v[1]
      if c < t:
        y = yMin
        x = p[0] + c * v[0]
        t = c
    elif v[1] > 0.0: # bottom
      if p[1] >= yMax: return nil
      c = (yMax - p[1]) / v[1]
      if c < t:
        y = yMax
        x = p[0] + c * v[0]
        t = c

    if v[0] > 0.0: # right
      if p[0] >= xMax: return nil
      c = (xMax - p[0]) / v[0]
      if c < t:
        x = xMax
        y = p[1] + c * v[1]
        t = c
    elif v[0] < 0.0: # left
      if p[0] <= xMin: return nil
      c = (xMin - p[0]) / v[0]
      if c < t:
        x = xMin
        y = p[1] + c * v[1]
        t = c

    return [x, y]

  #[
  function clockwise([x0, y0], [x1, y1], [x2, y2]) {
    return (x1 - x0) * (y2 - y0) < (y1 - y0) * (x2 - x0);
  }
  ]#
  func clockwise(p, q, r: array[2, T]): bool =
    return (q[0] - p[0]) * (r[1] - p[1]) < (q[1] - p[1]) * (r[0] - p[0])

  #[
  function contains({points, v0, vn}, p) {
    let n = points.length, p0, p1 = points[0];
    if (clockwise(p, [p1[0] + v0[0], p1[1] + v0[1]], p1)) return false;
    for (let i = 1; i < n; ++i) if (clockwise(p, p0 = p1, p1 = points[i])) return false;
    if (clockwise(p, p1, [p1[0] + vn[0], p1[1] + vn[1]])) return false;
    return true;
  }
  ]#
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

    #[
  inside = ({
    top: ([x, y]) => y > ymin,
    right: ([x, y]) => x < xmax,
    bottom: ([x, y]) => y < ymax,
    left: ([x, y]) => x > xmin
  })
  intersect = ({
    top: ([x0, y0], [x1, y1]) => [x0 + (x1 - x0) * (ymin - y0) / (y1 - y0), ymin],
    right: ([x0, y0], [x1, y1]) => [xmax, y0 + (y1 - y0) * (xmax - x0) / (x1 - x0)],
    bottom: ([x0, y0], [x1, y1]) => [x0 + (x1 - x0) * (ymax - y0) / (y1 - y0), ymax],
    left: ([x0, y0], [x1, y1]) => [xmin, y0 + (y1 - y0) * (xmin - x0) / (x1 - x0)]
  })
  ]#
  let
    insideTop = proc (p: array[2, T]): bool = p[1] > yMin
    insideRight = proc (p: array[2, T]): bool = p[0] < xMax
    insideBottom = proc (p: array[2, T]): bool = p[1] < yMax
    insideLeft = proc (p: array[2, T]): bool = p[0] > xMin

    intersectTop = proc (p, q: array[2, T]): array[2, T] = [p[0] + ( q[0] - p[0] ) * ( yMin - p[1] ) / ( q[1] - p[1] ), yMin]
    intersectRight = proc (p, q: array[2, T]): array[2, T] = [xMax, p[1] + ( q[1] - p[1] ) * ( xMax - p[0] ) / ( q[0] - p[0] )]
    intersectBottom = proc (p, q: array[2, T]): array[2, T] = [p[0] + ( q[0] - p[0] ) * ( yMax - p[1] ) / ( q[1] - p[1] ), yMax]
    intersectLeft = proc (p, q: array[2, T]): array[2, T] = [xMin, p[1] + ( q[1] - p[1] ) * ( xMin - p[0] ) / ( q[0] - p[0] )]


  #[
  function clipper(inside, intersect) {
    return function(subject) {
      const P = [], n = subject.length;
      if (!n) return P;
      let p0, p1 = subject[n - 1];
      let t0, t1 = inside(p1);
      for (let i = 0; i < n; ++i) {
        p0 = p1, p1 = subject[i];
        t0 = t1, t1 = inside(p1);
        if (t1 !== t0) P.push(intersect(p0, p1));
        if (t1) P.push(p1);
      }
      return P;
    };
  }
  ]#
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


  #[
  clip = {
    const top = clipper(inside.top, intersect.top);
    const right = clipper(inside.right, intersect.right);
    const bottom = clipper(inside.bottom, intersect.bottom);
    const left = clipper(inside.left, intersect.left);
    return function clip(subject) {
      return left(bottom(right(top(subject))));
    };
  }
  ]#
  let
    top = clipper(insideTop, intersectTop)
    right = clipper(insideRight, intersectRight)
    bottom = clipper(insideBottom, intersectBottom)
    left = clipper(insideLeft, intersectLeft)

  proc clip(subject: seq[array[2, T]]): seq[array[2, T]] =
    return left(bottom(right(top(subject))))

  #[
  function sidecode([x, y]) {
    return (x === xmin ? 1 : x === xmax ? 2 : 0) | (y === ymin ? 4 : y === ymax ? 8 : 0);
  }
  ]#
  func sidecode(p: array[2, T]): int8 =
    return (if p[0] == xMin: 1 else: (if p[0] == xMax: 2 else: 0)) or (if p[1] == yMin: 4 else: (if p[1] == yMax: 8 else: 0))

#[
  let P = polygon.points.slice(), p, n;
  if (p = project(P[0], polygon.v0)) P.unshift(p);
  if (p = project(P[P.length - 1], polygon.vn)) P.unshift(p);
  if (n = (P = clip(P)).length) {
    for (let i = 0, c0, c1 = sidecode(P[n - 1]); i < n; ++i) {
      c0 = c1, c1 = sidecode(P[i]);
      if (c0 && c1) {
        while (c0 !== c1) {
          let c;
          switch (c0) {
            case 0b0101: c0 = 0b0100; continue; // top-left
            case 0b0100: c0 = 0b0110, c = [xmax, ymin]; break; // top
            case 0b0110: c0 = 0b0010; continue; // top-right
            case 0b0010: c0 = 0b1010, c = [xmax, ymax]; break; // right
            case 0b1010: c0 = 0b1000; continue; // bottom-right
            case 0b1000: c0 = 0b1001, c = [xmin, ymax]; break; // bottom
            case 0b1001: c0 = 0b0001; continue; // bottom-left
            case 0b0001: c0 = 0b0101, c = [xmin, ymin]; break; // left
          }
          if (contains(polygon, c)) {
            P.splice(i, 0, c), ++n, ++i;
          }
        }
      }
    }
  } else if (contains(polygon, [(xmin + xmax) / 2, (ymin + ymax) / 2])) {
    P.push([xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax]);
  }
  return P;
]#
  var
    ply = polygon.points
    p: array[2, T]
    n: int
  p = project(ply[0], polygon.v0)
  if p: ply.insert(@[p]) # FIXME: Better way than insert?
  p = project(ply[^1], polygon.vn)
  if p: ply.insert(@[p])
  ply = clip(ply)
  n = ply.len
  if n > 0:
    var
      i = 0
      c0, c1: int8
    c1 = sidecode(ply[^1])
    while i < n:
      inc i
      c0 = c1
      c1 = sidecode(ply[i])
      if (c0 and c1) != 0:
        while c0 != c1:
          var c: array[2, T]
          case c0
          of 0b0101: # top-left
            c0 = 0b0100
            continue
          of 0b0100: # top
            c0 = 0b0110
            c = [xMax, yMin]
            break
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
            break
          of 0b1001: # bottom-left
            c0 = 0b0001
            continue
          of 0b0001: # left
            c0 = 0b0101
            c = [xMin, yMin]
            break
          if contains(polygon, c):
            ply.insert(@[c], i)
            inc n
            inc i

  elif contains(polygon, [(xMin + xMax) / 2.0, (yMin + yMax) / 2.0]):
    ply.insert(@[[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax]])

  return ply
