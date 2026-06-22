# lean-humanoid-rom

Humanoid range-of-motion / IK-constraint hexagon: ROM math, Kusudama, muscle/prismatic constraints (core); B3D/AddBiomechanics parsers + GPU shader + python tools (adapters). See `CITATIONS.bib` for the body-model / biomechanics / IK references.

> Split out of the [`lean-predictive-bvh`](https://github.com/v-sekai-multiplayer-fabric/lean-predictive-bvh) monorepo (now archived). Each hexagon cluster is its own repo following the `core/ports/adapters` convention; cross-cluster wiring is via Lake `require ... from git`.

## Dependencies

- [`lean-shared-core`](v-sekai-multiplayer-fabric/lean-shared-core) — common primitive types
- [`LeanSlang`](https://github.com/V-Sekai-fire/lean-slang) — Slang/HLSL AST for the Kusudama shader (pure-Lean use; no FFI link)

## Build

```sh
lake build         # production gate: typecheck the  cluster
lake build Research  # research-tier (non-gating; may fail)
```

## Hexagon layout

- `core/` — dependency-free domain logic + proofs
- `ports/` — narrow driving (source) / driven (sink) contracts
- `adapters/` — concrete I/O at the edges
