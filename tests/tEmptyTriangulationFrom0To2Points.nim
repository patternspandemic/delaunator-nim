
import ../src/delaunator

include "fixtures/ukraine.nim"

# No points
var
  coords: seq[float64] = @[]
  d = delaunator.fromCoords[float64](coords)
assert d.triangles.len == 0
assert d.hull.len == 0

# One points
d = delaunator.fromPoints[array[2, int], float64](ukraine[0..<1])
assert d.triangles.len == 0
assert d.hull == @[uint32(0)]

# Two points
d = delaunator.fromPoints[array[2, int], float64](ukraine[0..<2])
assert d.triangles.len == 0
assert d.hull == @[uint32(0), uint32(1)] or d.hull == @[uint32(1), uint32(0)]
