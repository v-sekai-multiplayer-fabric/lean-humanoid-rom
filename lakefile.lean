-- SPDX-License-Identifier: MIT
-- Copyright (c) 2026-present K. S. Ernest (iFire) Lee

import Lake
open System Lake DSL

package «lean-humanoid-rom» where

-- Shared primitive types (common vocabulary).
require «lean-shared-core» from git
  "https://github.com/v-sekai-multiplayer-fabric/lean-shared-core.git" @ "main"

-- LeanSlang: Lean 4 AST for Slang/HLSL shader emission (KusudamaShader adapter).
-- Used as a pure-Lean library here (no exe), so its libslang FFI is never linked
-- and the Slang SDK does not need vendoring for this repo's build.
require LeanSlang from git
  "https://github.com/V-Sekai-fire/lean-slang.git" @ "main"

-- Humanoid range-of-motion / IK-constraint hexagon: ROM math, Kusudama,
-- muscle/prismatic constraints (core); B3D/AddBiomechanics parsers + GPU shader
-- + python tools (adapters).
lean_lib HumanoidRom where
  roots := #[`HumanoidRom]
  globs := #[.one `HumanoidRom]

-- Research-tier ROMTool (NOT on the CI production gate; typeclass synth failure).
lean_lib Research where
  roots := #[`Research]
  globs := #[.one `Research]
