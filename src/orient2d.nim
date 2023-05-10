import std/math

let
  # from robust-predicates util.js
  epsilon = 1.1102230246251565e-16
  splitter = 134217729
  resulterrbound = (3 + 8 * epsilon) * epsilon
  # from robust-predicates orient2d.js
  ccwerrboundA = (3 + 16 * epsilon) * epsilon
  ccwerrboundB = (2 + 12 * epsilon) * epsilon
  ccwerrboundC = (9 + 64 * epsilon) * epsilon * epsilon

var
  # from robust-predicates orient2d.js
  # Direct array assignment instead of using vec func
  B: array[4, float64]
  C1: array[8, float64]
  C2: array[12, float64]
  D: array[16, float64]
  u: array[4, float64]


#[
func vec(n: int): seq[Float64] =
  return newSeq[Float64](n)
]#


func estimate(e: array[4, float64]): float64 = math.sum(e)


# fast_expansion_sum_zeroelim routine from oritinal code
proc sum(elen: int; e: openArray[float64]; flen: int; f: openArray[float64]; h: var openArray[float64]): int =
  var
    Q, Qnew, hh, bvirt: float64
    enow = e[0]
    fnow = f[0]
    eindex = 0
    findex = 0
  if (fnow > enow) == (fnow > -enow):
    Q = enow
    inc eindex
    enow = e[eindex]
  else:
    Q = fnow
    inc findex
    fnow = f[findex]
  var hindex = 0
  if eindex < elen and findex < flen:
    if (fnow > enow) == (fnow > -enow):
      Qnew = enow + Q
      hh = Q - (Qnew - enow)
      inc eindex
      enow = e[eindex]
    else:
      Qnew = fnow + Q
      hh = Q - (Qnew - fnow)
      inc findex
      fnow = f[findex]
    Q = Qnew
    if hh != 0:
      h[hindex] = hh
      inc hindex
    while eindex < elen and findex < flen:
      if (fnow > enow) == (fnow > -enow):
        Qnew = Q + enow
        bvirt = Qnew - Q
        hh = Q - (Qnew - bvirt) + (enow - bvirt)
        inc eindex
        enow = e[eindex]
      else:
        Qnew = Q + fnow
        bvirt = Qnew - Q
        hh = Q - (Qnew - bvirt) + (fnow - bvirt)
        inc findex
        fnow = f[findex]
      Q = Qnew
      if hh != 0:
        h[hindex] = hh
        inc hindex
  while eindex < elen:
    Qnew = Q + enow
    bvirt = Qnew - Q
    hh = Q - (Qnew - bvirt) + (enow - bvirt)
    inc eindex
    enow = e[eindex]
    Q = Qnew
    if hh != 0:
      h[hindex] = hh
      inc hindex
  while findex < flen:
    Qnew = Q + fnow
    bvirt = Qnew - Q
    hh = Q - (Qnew - bvirt) + (fnow - bvirt)
    inc findex
    fnow = f[findex]
    Q = Qnew
    if hh != 0:
      h[hindex] = hh
      inc hindex
  if Q != 0 or hindex == 0:
    h[hindex] = Q
    inc hindex
  return hindex


func orient2dadapt(ax, ay, bx, by, cx, cy, detsum: SomeFloat): int = discard


func orient2d*(ax, ay, bx, by, cx, cy: SomeFloat): int =
  let
    detleft = (ay - cy) * (bx - cx)
    detright = (ax - cx) * (by - cy)
    det = detleft - detright

  if detleft == 0 or detright == 0 or ((detleft > 0) != (detright > 0)):
    return det

  let detsum = abs(detleft + detright)
  if abs(det) >= ccwerrboundA * detsum:
    return det

  return -orient2dadapt(ax, ay, bx, by, cx, cy, detsum)
