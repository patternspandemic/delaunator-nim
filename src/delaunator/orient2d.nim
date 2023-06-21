# Up to date with mourner/robust-predicates c20b0ab9ab4c4f2969f3611908c41ce76aa0e7a7 May 25, 2023

import std/math

let
  # from robust-predicates util.js
  epsilon = 1.1102230246251565e-16
  splitter = 134217729.0
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


proc orient2dadapt(ax, ay, bx, by, cx, cy, detsum: SomeFloat): SomeFloat =
  var
    acxtail, acytail, bcxtail, bcytail: SomeFloat
    bvirt, c, ahi, alo, bhi, blo, o_i, o_j, o_0, s1, s0, t1, t0, u3: SomeFloat

  let
    acx = ax - cx
    bcx = bx - cx
    acy = ay - cy
    bcy = by - cy

  s1 = acx * bcy
  c = splitter * acx
  ahi = c - (c - acx)
  alo = acx - ahi
  c = splitter * bcy
  bhi = c - (c - bcy)
  blo = bcy - bhi
  s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo)
  t1 = acy * bcx
  c = splitter * acy
  ahi = c - (c - acy)
  alo = acy - ahi
  c = splitter * bcx
  bhi = c - (c - bcx)
  blo = bcx - bhi
  t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo)
  o_i = s0 - t0
  bvirt = s0 - o_i
  B[0] = s0 - (o_i + bvirt) + (bvirt - t0)
  o_j = s1 + o_i
  bvirt = o_j - s1
  o_0 = s1 - (o_j - bvirt) + (o_i - bvirt)
  o_i = o_0 - t1
  bvirt = o_0 - o_i
  B[1] = o_0 - (o_i + bvirt) + (bvirt - t1)
  u3 = o_j + o_i
  bvirt = u3 - o_j
  B[2] = o_j - (u3 - bvirt) + (o_i - bvirt)
  B[3] = u3

  var
    det = estimate(B)
    errbound = ccwerrboundB * detsum
  if det >= errbound or -det >= errbound:
     return det

  bvirt = ax - acx
  acxtail = ax - (acx + bvirt) + (bvirt - cx)
  bvirt = bx - bcx
  bcxtail = bx - (bcx + bvirt) + (bvirt - cx)
  bvirt = ay - acy
  acytail = ay - (acy + bvirt) + (bvirt - cy)
  bvirt = by - bcy
  bcytail = by - (bcy + bvirt) + (bvirt - cy)

  if acxtail == 0 and acytail == 0 and bcxtail == 0 and bcytail == 0:
    return det

  errbound = ccwerrboundC * detsum + resulterrbound * abs(det)
  det += (acx * bcytail + bcy * acxtail) - (acy * bcxtail + bcx * acytail)
  if det >= errbound or -det >= errbound:
    return det

  s1 = acxtail * bcy
  c = splitter * acxtail
  ahi = c - (c - acxtail)
  alo = acxtail - ahi
  c = splitter * bcy
  bhi = c - (c - bcy)
  blo = bcy - bhi
  s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo)
  t1 = acytail * bcx
  c = splitter * acytail
  ahi = c - (c - acytail)
  alo = acytail - ahi
  c = splitter * bcx
  bhi = c - (c - bcx)
  blo = bcx - bhi
  t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo)
  o_i = s0 - t0
  bvirt = s0 - o_i
  u[0] = s0 - (o_i + bvirt) + (bvirt - t0)
  o_j = s1 + o_i
  bvirt = o_j - s1
  o_0 = s1 - (o_j - bvirt) + (o_i - bvirt)
  o_i = o_0 - t1
  bvirt = o_0 - o_i
  u[1] = o_0 - (o_i + bvirt) + (bvirt - t1)
  u3 = o_j + o_i
  bvirt = u3 - o_j
  u[2] = o_j - (u3 - bvirt) + (o_i - bvirt)
  u[3] = u3
  let C1len = sum(4, B, 4, u, C1)

  s1 = acx * bcytail
  c = splitter * acx
  ahi = c - (c - acx)
  alo = acx - ahi
  c = splitter * bcytail
  bhi = c - (c - bcytail)
  blo = bcytail - bhi
  s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo)
  t1 = acy * bcxtail
  c = splitter * acy
  ahi = c - (c - acy)
  alo = acy - ahi
  c = splitter * bcxtail
  bhi = c - (c - bcxtail)
  blo = bcxtail - bhi
  t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo)
  o_i = s0 - t0
  bvirt = s0 - o_i
  u[0] = s0 - (o_i + bvirt) + (bvirt - t0)
  o_j = s1 + o_i
  bvirt = o_j - s1
  o_0 = s1 - (o_j - bvirt) + (o_i - bvirt)
  o_i = o_0 - t1
  bvirt = o_0 - o_i
  u[1] = o_0 - (o_i + bvirt) + (bvirt - t1)
  u3 = o_j + o_i
  bvirt = u3 - o_j
  u[2] = o_j - (u3 - bvirt) + (o_i - bvirt)
  u[3] = u3
  let C2len = sum(C1len, C1, 4, u, C2)

  s1 = acxtail * bcytail
  c = splitter * acxtail
  ahi = c - (c - acxtail)
  alo = acxtail - ahi
  c = splitter * bcytail
  bhi = c - (c - bcytail)
  blo = bcytail - bhi
  s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo)
  t1 = acytail * bcxtail
  c = splitter * acytail
  ahi = c - (c - acytail)
  alo = acytail - ahi
  c = splitter * bcxtail
  bhi = c - (c - bcxtail)
  blo = bcxtail - bhi
  t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo)
  o_i = s0 - t0
  bvirt = s0 - o_i
  u[0] = s0 - (o_i + bvirt) + (bvirt - t0)
  o_j = s1 + o_i
  bvirt = o_j - s1
  o_0 = s1 - (o_j - bvirt) + (o_i - bvirt)
  o_i = o_0 - t1
  bvirt = o_0 - o_i
  u[1] = o_0 - (o_i + bvirt) + (bvirt - t1)
  u3 = o_j + o_i
  bvirt = u3 - o_j
  u[2] = o_j - (u3 - bvirt) + (o_i - bvirt)
  u[3] = u3
  let Dlen = sum(C2len, C2, 4, u, D)

  return D[Dlen - 1]


proc orient2d*(ax, ay, bx, by, cx, cy: SomeFloat): SomeFloat =
  let
    detleft = (ay - cy) * (bx - cx)
    detright = (ax - cx) * (by - cy)
    det = detleft - detright

  #[ https://github.com/mourner/robust-predicates/pull/7
  # Remove redundant test for better branch predictability
  if detleft == 0 or detright == 0 or ((detleft > 0) != (detright > 0)):
    return det
  ]#

  let detsum = abs(detleft + detright)
  if abs(det) >= ccwerrboundA * detsum:
    return det

  return -orient2dadapt(ax, ay, bx, by, cx, cy, detsum)
