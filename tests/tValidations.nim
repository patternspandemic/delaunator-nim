
import std/sequtils

import ../src/delaunator
from ./helpers import validate

include "fixtures/ukraine.nim"
include "fixtures/robustness1.nim"
include "fixtures/robustness2.nim"
include "fixtures/robustness3.nim"
include "fixtures/robustness4.nim"
include "fixtures/issue13.nim"
include "fixtures/issue43.nim"
include "fixtures/issue44.nim"


validate[array[2, int], float64](ukraine)
validate[array[2, int], float64](@[[516, 661], [369, 793], [426, 539], [273, 525], [204, 694], [747, 750], [454, 390]]) # issue11
validate[array[2, float64], float64](issue13)
validate[array[2, int], float64](@[[382, 302], [382, 328], [382, 205], [623, 175], [382, 188], [382, 284], [623, 87], [623, 341], [141, 227]]) # issue 24
validate[array[2, float64], float64](issue43)
validate[array[2, float64], float64](issue44)
validate[array[2, float64], float64](robustness1)
validate[array[2, float64], float64](robustness1.map(proc (p: array[2, float64]): array[2, float64] = [p[0] / 1e9, p[1] / 1e9]))
validate[array[2, float64], float64](robustness1.map(proc (p: array[2, float64]): array[2, float64] = [p[0] / 100.0, p[1] / 100.0]))
validate[array[2, float64], float64](robustness1.map(proc (p: array[2, float64]): array[2, float64] = [p[0] * 100.0, p[1] * 100.0]))
validate[array[2, float64], float64](robustness1.map(proc (p: array[2, float64]): array[2, float64] = [p[0] * 1e9, p[1] * 1e9]))
validate[array[2, float64], float64](robustness2[0 ..< 100])
validate[array[2, float64], float64](robustness2)
validate[array[2, float64], float64](robustness3)
validate[array[2, float64], float64](robustness4)
