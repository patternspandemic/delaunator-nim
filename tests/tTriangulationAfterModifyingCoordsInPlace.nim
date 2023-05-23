
import std/sequtils

import ../src/delaunator
from ./helpers import validate

include "fixtures/ukraine.nim"


var d = delaunator.fromPoints[array[2, int], float64](ukraine)

validate[array[2, int], float64](ukraine, d)
assert d.triangles.len == 5133

let p: array[2, int] = [80, 220]
d.coords[0] = float64(p[0])
d.coords[1] = float64(p[1])

var newPoints: seq[array[2, int]] = concat[array[2, int]](@[p], ukraine[1..^1])

update[float64](d)
validate[array[2, int], float64](newPoints, d)
assert d.triangles.len == 5139
