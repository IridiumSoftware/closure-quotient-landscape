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
  , q181SeedEdges
  , q181NVerts
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

-- | The canonical Q_181 seed: 132-edge irregular topology on 12 vertices,
-- baked here as a literal so |Q_181| = 181 is reproducible without re-running
-- the random density sweep. Discovered by 'paper/scripts/q181_search.py'
-- under @random.seed(2718)@, density 0.10, trial 0; the discovery method is
-- the search but the existence claim rests on this fixed edge list.
q181NVerts :: Int
q181NVerts = 12

q181SeedEdges :: [Hyperedge]
q181SeedEdges =
  [ (0,2,1), (0,4,9), (0,4,10), (0,8,3), (0,9,1), (0,9,7)
  , (0,9,8), (0,10,4), (0,11,5), (1,0,9), (1,2,10), (1,4,7)
  , (1,6,4), (1,6,8), (1,7,4), (1,7,5), (1,8,0), (1,10,4)
  , (2,0,6), (2,1,9), (2,3,1), (2,3,5), (2,3,6), (2,6,1)
  , (2,6,10), (2,7,1), (2,9,1), (2,10,5), (3,0,10), (3,0,11)
  , (3,2,10), (3,4,5), (3,4,7), (3,5,6), (3,7,5), (3,10,1)
  , (4,0,5), (4,0,8), (4,1,5), (4,1,10), (4,1,11), (4,2,9)
  , (4,3,1), (4,5,2), (4,5,9), (4,7,8), (4,10,11), (4,11,3)
  , (4,11,7), (5,1,4), (5,2,6), (5,2,7), (5,3,4), (5,4,10)
  , (5,6,4), (5,8,10), (5,8,11), (5,9,6), (5,9,8), (5,9,11)
  , (5,11,0), (6,0,5), (6,0,8), (6,0,9), (6,2,0), (6,2,7)
  , (6,4,9), (6,5,1), (6,7,0), (6,7,3), (6,7,4), (6,9,7)
  , (6,10,7), (7,0,3), (7,2,1), (7,2,9), (7,2,11), (7,3,4)
  , (7,5,10), (7,6,9), (7,6,10), (7,6,11), (7,10,0), (7,11,1)
  , (8,1,3), (8,1,5), (8,1,9), (8,2,3), (8,4,7), (8,4,9)
  , (8,5,4), (8,6,1), (8,6,7), (8,7,9), (8,9,10), (8,11,9)
  , (9,0,5), (9,0,8), (9,1,3), (9,2,5), (9,3,0), (9,4,10)
  , (9,5,7), (9,5,10), (9,5,11), (9,6,7), (9,7,10), (9,8,0)
  , (9,10,1), (9,11,7), (10,1,4), (10,1,8), (10,1,11), (10,2,7)
  , (10,2,9), (10,3,1), (10,3,6), (10,7,2), (10,8,0), (10,8,2)
  , (10,9,11), (10,11,1), (10,11,6), (11,1,8), (11,2,6), (11,3,9)
  , (11,4,2), (11,4,3), (11,5,2), (11,6,0), (11,9,10), (11,10,6)
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
