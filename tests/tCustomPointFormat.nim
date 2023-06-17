
import ../src/delaunator

type
  Site = tuple
    label: string
    x, y: int

let
  myGetX = proc (t: Site): float64 = float64(t.x)
  myGetY = proc (t: SIte): float64 = float64(t.y)
var
  sites: seq[Site] = @[("a", 5, 5), ("b", 7, 5), ("c", 7, 6)]
  dCustom = delaunator.fromCustom[Site, float64](sites, myGetX, myGetY)

assert dCustom.triangles == @[uint32(0), uint32(2), uint32(1)]

var
  points: seq[seq[int]] = @[@[5, 5], @[7, 5], @[7, 6]]
  dPoints = delaunator.fromPoints[seq[int], float32](points)

assert dPoints.triangles == @[uint32(0), uint32(2), uint32(1)]
