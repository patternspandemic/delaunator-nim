
import std/math
import ../src/delaunator


func orient[T](p, r, q: array[2, SomeNumber]): T =
  let
    l = T((r[1] - p[1]) * (q[0] - p[0]))
    r = T((r[0] - p[0]) * (q[1] - p[1]))
  if abs(l - r) >= 3.3306690738754716e-16 * abs(l + r):
    return l - r
  else:
    return 0.0


func convex[T](r, q, p: array[2,SomeNumber]): bool =
  # js equivelant of:
  # return (orient(p, r, q) || orient(r, q, p) || orient(q, p, r)) >= 0;
  var o = orient[T](p, r, q)
  if o != 0: return o >= 0
  o = orient[T](r, q, p)
  if o != 0: return o >= 0
  return orient[T](q, p, r) >= 0


# Kahan and Babuska summation, Neumaier variant; accumulates less FP error
func sum(x: seq[SomeFloat]): SomeFloat =
  var
    sum = x[0]
    err = 0.0
  for i in 1 ..< x.len:
    let
      k = x[i]
      m = sum + k
    let tmp = if abs(sum) >= abs(k): sum - m + k else: k - m + sum
    err += tmp
    sum = m
  return sum + err


proc validate*[P, T](points: seq[P], d: Delaunator[T]) =
  # validate halfedges
  for i in 0 ..< d.halfedges.len:
    assert (d.halfedges[i] == int32(-1) or d.halfedges[d.halfedges[i]] == int32(i)), "valid halfedge connection"

  # validate triangulation
  var
    hullAreas = newSeq[T]()
    hlen = d.hull.len
    j = hlen - 1
  for i in 0 ..< hlen:
    let
      tmp0 = points[d.hull[j]]
      tmp = points[d.hull[i]]
    hullAreas.add(T((tmp[0] - tmp0[0]) * (tmp[1] + tmp0[1])))
    assert convex[T](points[d.hull[j]], points[d.hull[(j + 1) mod hlen]],  points[d.hull[(j + 3) mod hlen]]), "hull should be convex at " & $j
    j = i
  let hullArea = sum(hullAreas)

  var
    triangleAreas = newSeq[T]()
    i = 0
  while i < d.triangles.len:
    let
      a = points[d.triangles[i]]
      b = points[d.triangles[i + 1]]
      c = points[d.triangles[i + 2]]
    triangleAreas.add(T(abs((b[1] - a[1]) * (c[0] - b[0]) - (b[0] - a[0]) * (c[1] - b[1]))))
    i += 3
  let trianglesArea = sum(triangleAreas)

  let err = abs((hullArea - trianglesArea) / hullArea)
  assert err <= pow(2.0, -51.0), "triangulation should be valid; " & $err & " error"


proc validate*[P, T](points: seq[P]) =
  validate(points, delaunator.fromPoints[P, T](points))
