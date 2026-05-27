-- |
-- Module      : CFS.Topologies
-- Description : Canonical seed topologies and Gaussian-integer initial conditions.
--
-- Mirrors the Python reference (@q102_build_exact_v1.py@): the same
-- topology builders (@G0_TOPO@, @adjacent_ternary@, @complete_ternary@)
-- and the same canonical 15-element @_GEN_ICS@ Gaussian-integer initial
-- conditions verified to put @K_6^3 -> Q_102@ exactly.
module CFS.Topologies
  ( -- * Seed topologies
    g0Topo
  , adjacentTernary
  , completeTernary
  , q45SeedEdges
    -- * Canonical Gaussian-integer ICs
  , canonicalICs
  ) where

import           CFS.Cochain
import           CFS.Multiway        (Hyperedge, Vertex)
import qualified Data.List           as List
import qualified Data.Map.Strict     as Map
import           Data.Map.Strict     (Map)

-- | The G_0 6-vertex Rosen-closed seed: |Q| = 24 single side, |Q| = 48
-- C-closed.
g0Topo :: [Hyperedge]
g0Topo =
  [(0, 1, 2), (1, 2, 3), (2, 3, 4), (3, 4, 5), (4, 5, 0), (5, 0, 1)]

-- | The adjacent-triple seed on 6 vertices: all permutations of
-- @(i, i+1 mod 6, i+2 mod 6)@. |seed| = 36; |Q| = 42 single side, |Q| = 84
-- C-closed.
adjacentTernary :: [Hyperedge]
adjacentTernary = List.sort $ List.nub
  [ p
  | i <- [0 .. 5 :: Int]
  , p <- permutations3 (i, (i + 1) `mod` 6, (i + 2) `mod` 6)
  ]

permutations3 :: Hyperedge -> [Hyperedge]
permutations3 (a, b, c) =
  [(a, b, c), (a, c, b), (b, a, c), (b, c, a), (c, a, b), (c, b, a)]

-- | The complete-ternary graph on @n@ vertices: all ordered triples of
-- distinct vertices. @K_6^3@ has |seed| = 120; |Q| = 51 single side,
-- |Q| = 102 C-closed.
completeTernary :: Int -> [Hyperedge]
completeTernary n =
  [ (i, j, k)
  | i <- [0 .. n - 1]
  , j <- [0 .. n - 1], j /= i
  , k <- [0 .. n - 1], k /= i, k /= j
  ]

-- | The Q_45 / Q_90 seed: 'adjacentTernary' plus six pads sampled from
-- @K_6^3 \\ adjacentTernary@ with @numpy.random.default_rng(42)@.
-- The pads are baked in here as a static list to keep the Haskell port
-- numpy-free; they reproduce the Python @q45_seed_edges(seed=42)@ exactly.
-- |Q| = 45 single side, |Q| = 90 C-closed.
q45SeedEdges :: [Hyperedge]
q45SeedEdges =
  adjacentTernary
    ++ [ (0, 2, 5)
       , (4, 2, 0)
       , (0, 2, 3)
       , (4, 1, 0)
       , (5, 3, 2)
       , (1, 5, 4)
       ]

-- | The canonical Gaussian-integer initial conditions verified to put
-- @K_6^3 -> Q_102@ exactly under the binary-operator / ternary-hyperedge
-- rewriting (mirrors the Python @_GEN_ICS@ list, first @n@ entries).
canonicalICs :: Int -> Map Vertex Cochain
canonicalICs n = Map.fromList (zip [0 ..] (take n genICs))

genICs :: [Cochain]
genICs =
  [ mkv ( 2,  1) ( 1,  0) ( 3, -1)
  , mkv ( 1,  0) ( 2,  1) ( 1, -2)
  , mkv ( 1, -1) ( 3,  0) ( 2,  1)
  , mkv ( 3,  0) ( 1, -1) ( 1,  2)
  , mkv ( 1,  2) ( 2, -1) ( 1,  0)
  , mkv ( 2,  0) ( 1,  1) ( 3,  0)
  , mkv ( 3,  1) ( 2,  0) ( 1, -1)
  , mkv ( 1,  1) ( 3, -1) ( 2,  0)
  , mkv ( 2, -1) ( 1,  2) ( 3,  1)
  , mkv ( 1,  2) ( 3,  1) ( 2, -1)
  , mkv ( 2, -1) ( 1,  3) ( 1,  1)
  , mkv ( 3,  1) ( 2, -1) ( 1,  2)
  , mkv ( 1,  1) ( 2,  2) ( 3, -1)
  , mkv ( 2,  1) ( 1, -2) ( 2,  1)
  , mkv ( 1, -1) ( 3,  2) ( 1,  1)
  ]
  where
    mkv x y z = (toG x, toG y, toG z)
    toG (r, i) = Gauss r i
