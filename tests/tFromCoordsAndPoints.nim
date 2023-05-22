
import ../src/delaunator

include "fixtures/normal50.nim"

var
  dFromCoords = delaunator.fromCoords[float64](normalCoords)
  dFromPoints = delaunator.fromPoints[array[2, float64], float64](normalPoints)

assert dFromCoords.coords == dFromPoints.coords
assert dFromCoords.triangles == dFromPoints.triangles
assert dFromCoords.halfedges == dFromPoints.halfedges
assert dFromCoords.hull == dFromPoints.hull
assert dFromCoords.minX == dFromPoints.minX
assert dFromCoords.minY == dFromPoints.minY
assert dFromCoords.maxX == dFromPoints.maxX
assert dFromCoords.maxY == dFromPoints.maxY

# TODO: Not sure how to test this. Using errormsg and file in testament spec didn't work.
# Triggers a {.error: "Coordinates must be seq[float32] or seq[float64] but got ..." pragma.
#discard delaunator.fromPoints[array[2, int], int](@[[1, 2],[3, 4],[5,6]])
