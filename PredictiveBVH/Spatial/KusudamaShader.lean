-- SPDX-License-Identifier: MIT
-- Copyright (c) 2026-present K. S. Ernest (iFire) Lee
--
-- Discrete Lipschitz verification for kusudama polygon projection.
-- Discretize S² via the 6 octahedron vertices + 12 edge midpoints (18 pts).
-- Check every adjacent pair: output distSq ≤ input distSq (non-expansive).

import PredictiveBVH.Spatial.SphericalPolygon

namespace PredictiveBVH.KusudamaShader

open PredictiveBVH.SphericalPolygon

private def sc : Int := 100

-- ── Octahedral vertices: 6 axis-aligned + 12 edge midpoints = 18 points ────

private def octaPoints : Array Vec3 := #[
  -- 6 axis vertices
  { x :=  sc, y :=  0,  z :=  0 },   -- 0: +X
  { x := -sc, y :=  0,  z :=  0 },   -- 1: -X
  { x :=  0,  y :=  sc, z :=  0 },   -- 2: +Y
  { x :=  0,  y := -sc, z :=  0 },   -- 3: -Y
  { x :=  0,  y :=  0,  z :=  sc },  -- 4: +Z
  { x :=  0,  y :=  0,  z := -sc },  -- 5: -Z
  -- 12 edge midpoints (normalized to sc via ≈ 70 ≈ sc/√2)
  { x :=  70, y :=  70, z :=  0 },   -- 6:  +X+Y
  { x :=  70, y := -70, z :=  0 },   -- 7:  +X-Y
  { x := -70, y :=  70, z :=  0 },   -- 8:  -X+Y
  { x := -70, y := -70, z :=  0 },   -- 9:  -X-Y
  { x :=  70, y :=  0,  z :=  70 },  -- 10: +X+Z
  { x :=  70, y :=  0,  z := -70 },  -- 11: +X-Z
  { x := -70, y :=  0,  z :=  70 },  -- 12: -X+Z
  { x := -70, y :=  0,  z := -70 },  -- 13: -X-Z
  { x :=  0,  y :=  70, z :=  70 },  -- 14: +Y+Z
  { x :=  0,  y :=  70, z := -70 },  -- 15: +Y-Z
  { x :=  0,  y := -70, z :=  70 },  -- 16: -Y+Z
  { x :=  0,  y := -70, z := -70 }   -- 17: -Y-Z
]

-- Adjacency: octahedral edges (vertex to midpoint, midpoint to vertex).
-- Each axis vertex connects to 4 midpoints; each midpoint connects to 2 axis verts.
private def octaEdges : Array (Nat × Nat) := #[
  -- +X edges
  (0, 6), (0, 7), (0, 10), (0, 11),
  -- -X edges
  (1, 8), (1, 9), (1, 12), (1, 13),
  -- +Y edges
  (2, 6), (2, 8), (2, 14), (2, 15),
  -- -Y edges
  (3, 7), (3, 9), (3, 16), (3, 17),
  -- +Z edges
  (4, 10), (4, 12), (4, 14), (4, 16),
  -- -Z edges
  (5, 11), (5, 13), (5, 15), (5, 17),
  -- Face-interior edges (midpoint to midpoint on same face)
  (6, 10), (6, 14), (7, 10), (7, 16),
  (8, 12), (8, 14), (9, 12), (9, 16),
  (6, 11), (6, 15), (7, 11), (7, 17),
  (8, 13), (8, 15), (9, 13), (9, 17),
  (10, 14), (10, 16), (12, 14), (12, 16),
  (11, 15), (11, 17), (13, 15), (13, 17)
]

-- ── Solver + violation counting ─────────────────────────────────────────────

/-- Angular distance proxy on unnormalized integer vectors.
    Uses `1 - cos(angle) = 1 - (a·b)/(|a||b|)` scaled to avoid division:
    `|a|²|b|² - (a·b)²` is zero when parallel, positive when diverging.
    This is monotone with angle² for angles in [0, π] and works on
    unnormalized integer vectors without any division or sqrt. -/
private def angDistSq (a b : Vec3) : Int :=
  let ab := dot a b
  let aa := dot a a
  let bb := dot b b
  aa * bb - ab * ab

private def solve (p : Vec3) (poly : ConvexPolygon) : Vec3 :=
  if insidePoly p poly then p else projectToPoly p poly

/-- Scale-invariant non-expansive check.
    sin²(angle_out) ≤ sin²(angle_in) iff
    angDistSq(out) * |in_i|² * |in_j|² ≤ angDistSq(in) * |out_i|² * |out_j|²
    where angDistSq(a,b) = |a|²|b|² - (a·b)² = |a|²|b|²sin²(θ). -/
private def isExpansive (pi pj oi oj : Vec3) : Bool :=
  let inAng  := angDistSq pi pj
  let outAng := angDistSq oi oj
  let inNormSq  := dot pi pi * dot pj pj
  let outNormSq := dot oi oi * dot oj oj
  -- outAng / outNormSq > inAng / inNormSq  ⟺  outAng * inNormSq > inAng * outNormSq
  if outNormSq == 0 || inNormSq == 0 then false
  else outAng * inNormSq > inAng * outNormSq

private def countViolations (poly : ConvexPolygon) : Nat :=
  let outs := octaPoints.map fun p => solve p poly
  octaEdges.foldl (fun acc (i, j) =>
    if isExpansive octaPoints[i]! octaPoints[j]! outs[i]! outs[j]!
    then acc + 1 else acc
  ) 0

-- ── Test polygons ───────────────────────────────────────────────────────────

private def squarePoly : ConvexPolygon :=
  mkPolygon #[
    { x :=  sc, y :=  0,  z := 0 },
    { x :=  0,  y :=  sc, z := 0 },
    { x := -sc, y :=  0,  z := 0 },
    { x :=  0,  y := -sc, z := 0 }
  ] { x := 0, y := 0, z := sc }

private def trianglePoly : ConvexPolygon :=
  mkPolygon #[
    { x :=  sc, y :=  0,  z := 0 },
    { x :=  0,  y :=  sc, z := 0 },
    { x :=  0,  y :=  0,  z := sc }
  ] { x := sc, y := sc, z := sc }

private def scrambledPoly : ConvexPolygon :=
  mkPolygon #[
    { x :=  sc, y :=  0,  z := 0 },
    { x := -sc, y :=  0,  z := 0 },
    { x :=  0,  y :=  sc, z := 0 },
    { x :=  0,  y := -sc, z := 0 }
  ] { x := 0, y := 0, z := sc }

-- ── Find violations ─────────────────────────────────────────────────────────

private def reportViolations (name : String) (poly : ConvexPolygon) : IO Unit := do
  let outs := octaPoints.map fun p => solve p poly
  let mut violations := 0
  for h : idx in [:octaEdges.size] do
    let (i, j) := octaEdges[idx]
    if isExpansive octaPoints[i]! octaPoints[j]! outs[i]! outs[j]! then
      violations := violations + 1
      IO.println s!"{name} TELEPORT edge ({i},{j})"
      IO.println s!"  in:  ({octaPoints[i]!.x},{octaPoints[i]!.y},{octaPoints[i]!.z}) -> ({octaPoints[j]!.x},{octaPoints[j]!.y},{octaPoints[j]!.z})"
      IO.println s!"  out: ({outs[i]!.x},{outs[i]!.y},{outs[i]!.z}) -> ({outs[j]!.x},{outs[j]!.y},{outs[j]!.z})"
  IO.println s!"{name}: {violations} violations"

#eval! do
  reportViolations "square" squarePoly
  reportViolations "triangle" trianglePoly
  reportViolations "scrambled" scrambledPoly

-- ── Coverage check ──────────────────────────────────────────────────────────

theorem point_count : octaPoints.size = 18 := by native_decide
theorem edge_count  : octaEdges.size = 48  := by native_decide

end PredictiveBVH.KusudamaShader
