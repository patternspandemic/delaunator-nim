
import ../src/delaunator

var
  dFromCoords = delaunator.fromCoords[float64](@[])
  dFromPoints = delaunator.fromPoints[array[2, float64], float64](@[])

assert dFromCoords.coords.len == 0
assert dFromCoords.triangles.len == 0
assert dFromCoords.halfedges.len == 0
assert dFromCoords.hull.len == 0
assert dFromCoords.minX == Inf
assert dFromCoords.minY == Inf
assert dFromCoords.maxX == NegInf
assert dFromCoords.maxY == NegInf
