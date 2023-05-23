
import ../src/delaunator
from ./helpers import validate

include "fixtures/issue13.nim"

validate[array[2, float64], float64](issue13)
