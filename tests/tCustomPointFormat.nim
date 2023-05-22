
import ../src/delaunator

var
  points: seq[tuple[x,y: int]] = @[(x: 5, y: 5), (x: 7, y: 5), (x: 7, y: 6)]
  d = delaunator.fromPoints[tuple[x, y: int], float64](points, proc (p: tuple[x,y: int]): float64 = float64(p.x), proc (p: tuple[x,y: int]): float64 = float64(p.y))

assert d.triangles == @[uint32(0), uint32(2), uint32(1)]
