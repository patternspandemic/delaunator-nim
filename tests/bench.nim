
import std/[math, monotimes, random, sequtils, strutils, terminal, times]
import ../src/delaunator


proc uniform[T](count: int): seq[array[2, T]] =
  var points: seq[array[2, T]] = newSeqOfCap[array[2, T]](count)
  for i in 0 ..< count:
    points.add([rand[T](T(0.0) .. T(1.0)) * T(1e3), rand[T](T(0.0) .. T(1.0)) * T(1e3)])
  return points


proc gaussian[T](count: int): seq[array[2, T]] =
  var points: seq[array[2, T]] = newSeqOfCap[array[2, T]](count)
  for i in 0 ..< count:
    points.add([T(gauss() * 1e3), T(gauss() * 1e3)])
  return points


proc grid[T](count: int): seq[array[2, T]] =
  var points: seq[array[2, T]] = newSeqOfCap[array[2, T]](count)
  let size = sqrt(T(count))
  for i in 0 ..< int(floor(size)):
    for j in 0 ..< int(floor(size)):
      points.add([T(i), T(j)])
  return points


proc degenerate[T](count: int): seq[array[2, T]] =
  var points: seq[array[2, T]] = newSeqOfCap[array[2, T]](count + 1)
  points.add([T(0.0), T(0.0)])
  for i in 0 ..< count:
    let angle = 2.0 * PI * T(i) / T(count)
    points.add([T(1e10 * sin(angle)), T(1e10 * cos(angle))])
  return points


proc triangulateCoords[T](points: seq[T]): Delaunator[T] =
  return delaunator.fromCoords[T](points)


proc triangulatePoints[P, T](points: seq[array[2, T]]): Delaunator[T] =
  return delaunator.fromPoints[P, T](points)


when isMainModule:
  let
    distNames = @["uniform", "gaussian", "grid", "degenerate"]
    f32Distributions = @[uniform[float32], gaussian[float32], grid[float32], degenerate[float32]]
    f64Distributions = @[uniform[float64], gaussian[float64], grid[float64], degenerate[float64]]
    counts = @[20_000, 100_000, 200_000, 500_000, 1_000_000]

  randomize()

  echo ""
  writeStyled("      ")
  writeStyled("n-points:", {styleDim, styleUnderscore, styleItalic})
  writeStyled("   20k    100k    200k    500k   1000k\n", {styleDim, styleUnderscore, styleItalic})
  echo ""

  # bench float32 points
  writeStyled("  fromPoints (float32):\n", {styleBright})
  for i in 0 ..< distNames.len:
    writeStyled("    " & align(distNames[i], 10) & ":", {styleDim})

    let generate = f32Distributions[i]

    # warmup
    discard triangulatePoints[array[2, float32], float32](generate(counts[0]))
    discard triangulatePoints[array[2, float32], float32](generate(counts[1]))

    for c in counts:
      var points = generate(c)
      let strt = getMonotime()
      discard triangulatePoints[array[2, float32], float32](points)
      let elpsd = (getMonotime() - strt).inMilliseconds
      writeStyled(align($elpsd & "ms", 6) & "  ", {styleDim})
      flushFile(stdout)
    writeLine(stdout, "")

  # bench float64 points
  writeLine(stdout, "")
  writeStyled("  fromPoints (float64):\n", {styleBright})
  for i in 0 ..< distNames.len:
    writeStyled("    " & align(distNames[i], 10) & ":", {styleDim})

    let generate = f64Distributions[i]

    # warmup
    discard triangulatePoints[array[2, float64], float64](generate(counts[0]))
    discard triangulatePoints[array[2, float64], float64](generate(counts[1]))

    for c in counts:
      var points = generate(c)
      let strt = getMonotime()
      discard triangulatePoints[array[2, float64], float64](points)
      let elpsd = (getMonotime() - strt).inMilliseconds
      writeStyled(align($elpsd & "ms", 6) & "  ", {styleDim})
      flushFile(stdout)
    writeLine(stdout, "")

  # bench float32 coords
  writeLine(stdout, "")
  writeStyled("  fromCoords (float32):\n", {styleBright})
  for i in 0 ..< distNames.len:
    writeStyled("    " & align(distNames[i], 10) & ":", {styleDim})

    let generate = f32Distributions[i]

    var
      pointsC0 = generate(counts[0])
      coordsC0 = newSeq[float32](pointsC0.len * 2)
      pointsC1 = generate(counts[1])
      coordsC1 = newSeq[float32](pointsC1.len * 2)
    for i, point in pointsC0:
      coordsC0[2 * i] = point[0]
      coordsC0[2 * i + 1] = point[1]
    for i, point in pointsC1:
      coordsC1[2 * i] = point[0]
      coordsC1[2 * i + 1] = point[1]

    # warmup
    discard triangulateCoords[float32](coordsC0)
    discard triangulateCoords[float32](coordsC1)

    for c in counts:
      var
        points = generate(c)
        coords = newSeq[float32](points.len * 2)
      for i, point in points:
        coords[2 * i] = point[0]
        coords[2 * i + 1] = point[1]

      let strt = getMonotime()
      discard triangulateCoords[float32](coords)
      let elpsd = (getMonotime() - strt).inMilliseconds
      writeStyled(align($elpsd & "ms", 6) & "  ", {styleDim})
      flushFile(stdout)
    writeLine(stdout, "")

  # bench float64 coords
  writeLine(stdout, "")
  writeStyled("  fromCoords (float64):\n", {styleBright})
  for i in 0 ..< distNames.len:
    writeStyled("    " & align(distNames[i], 10) & ":", {styleDim})

    let generate = f64Distributions[i]

    var
      pointsC0 = generate(counts[0])
      coordsC0 = newSeq[float64](pointsC0.len * 2)
      pointsC1 = generate(counts[1])
      coordsC1 = newSeq[float64](pointsC1.len * 2)
    for i, point in pointsC0:
      coordsC0[2 * i] = point[0]
      coordsC0[2 * i + 1] = point[1]
    for i, point in pointsC1:
      coordsC1[2 * i] = point[0]
      coordsC1[2 * i + 1] = point[1]

    # warmup
    discard triangulateCoords[float64](coordsC0)
    discard triangulateCoords[float64](coordsC1)

    for c in counts:
      var
        points = generate(c)
        coords = newSeq[float64](points.len * 2)
      for i, point in points:
        coords[2 * i] = point[0]
        coords[2 * i + 1] = point[1]

      let strt = getMonotime()
      discard triangulateCoords[float64](coords)
      let elpsd = (getMonotime() - strt).inMilliseconds
      writeStyled(align($elpsd & "ms", 6) & "  ", {styleDim})
      flushFile(stdout)
    writeLine(stdout, "")

  stdout.resetAttributes()
