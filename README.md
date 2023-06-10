# Delaunator
Fast 2D Delaunay triangulation. A Nim port of [Mapbox/Delaunator](https://github.com/mapbox/delaunator).

> Note: This port of Delaunator is not optimized. That said, give the benchmark a run an see if it measures up. Even better, help optimize the code! ([notes on this port](#notes-on-this-port))

`nimble install delaunator`

[API Reference](https://patternspandemic.github.io/delaunator-nim/) - A work in progress.

See also, [Pixienator](https://github.com/patternspandemic/pixienator), a helper library for visualizing Delaunator using [Pixie](https://github.com/treeform/pixie).

<img src="images/delaunator.png" alt="Delaunator generated image example.">

### Features
- Delaunay Triangulation
- Voronoi Regions
- Uses porteded robust-predicates for robustness
- Helpers for navigating various parts of the datastructure (not optimized)
  - Includes clipping of infinite regions, ala d3-delaunay's implementation

### Example
```nim
# Example coming soon..
```

### Performance
I'd post some numbers here, but my kit is so old, you best just run the benchmark yourself. See [tests/bench.nim](https://github.com/patternspandemic/delaunator-nim/blob/main/tests/bench.nim) for a benchmark based off the one in the original implementation.

### Notes on This Port
This port has been implemented with a novice understanding of Nim, with all that entails. Practically speaking, this means the code has room for optimization, may not adhere to nim conventions and styling, and may contain outright puzzling verbosity / actions.

Some of this is on purpose. For instance, keeping code structure and variable naming as close to the original implementation made it easier to perform the port and locate bugs.

I've attempted to support both floating point types for coordinates via generics. This appears to have made the code a bit rough around the edges with casting here and there to appease the compiler.

Furthermore, the way halfedges are used as indices into triangles (and therefor coordinates) is problematic (I think) to access the full range of points beyond the high(int32) without changing halfedges to be int64. More testing will have to be done to determine tradeoffs.

In general, I suspect a number of edge cases have yet to be handled.
