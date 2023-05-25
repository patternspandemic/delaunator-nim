
import ../src/delaunator

var
  coords: seq[float64] = @[0.0, 0.0, 1.0, 0.0, 3.0, 0.0, 2.0, 0.0]
  d = delaunator.fromCoords[float64](coords)
assert d.triangles.len == 0
assert d.hull == @[uint32(0), uint32(1), uint32(3), uint32(2)] or d.hull == @[uint32(2), uint32(3), uint32(0), uint32(1)]
