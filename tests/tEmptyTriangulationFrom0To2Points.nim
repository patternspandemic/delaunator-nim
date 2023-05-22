
import ../src/delaunator

include "fixtures/ukraine.nim"

# No points
var d = delaunator.fromCoords[float64](@[])
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
