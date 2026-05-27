{-# LANGUAGE BangPatterns      #-}
{-# LANGUAGE OverloadedStrings #-}
-- |
-- Module      : CFS.Cochain
-- Description : Gaussian-integer 3-vectors and the binary conjugated cross-product.
--
-- The CFS rewriting rule's /operator/ is binary:
--
-- > compose :: Cochain -> Cochain -> Cochain
--
-- That is the type-level statement of the operator's arity. The ternary
-- structure lives one level up, in 'CFS.Multiway.Hyperedge', not here.
module CFS.Cochain
  ( -- * Gaussian integers
    Gauss(..)
  , gmul, gsub
  , gconj
  , gisZero
    -- * Cochains: 3-vectors over Z[i]
  , Cochain
  , cross3
  , compose
  , projEquiv
  , isZeroVec
  , conjVec
  ) where

-- | Gaussian integer (re :+ im).
data Gauss = Gauss !Integer !Integer
  deriving (Eq, Ord, Show)

gmul :: Gauss -> Gauss -> Gauss
gmul (Gauss a b) (Gauss c d) = Gauss (a*c - b*d) (a*d + b*c)

gsub :: Gauss -> Gauss -> Gauss
gsub (Gauss a b) (Gauss c d) = Gauss (a - c) (b - d)

gconj :: Gauss -> Gauss
gconj (Gauss a b) = Gauss a (-b)

gisZero :: Gauss -> Bool
gisZero (Gauss a b) = a == 0 && b == 0

-- | A cochain is a Gaussian-integer 3-vector. The triple is the
-- /element type/, not a relational tuple; the rewriting graph's ternary
-- hyperedges live in 'CFS.Multiway'.
type Cochain = (Gauss, Gauss, Gauss)

-- | Standard complex cross product, extended to Z[i].
cross3 :: Cochain -> Cochain -> Cochain
cross3 (a0, a1, a2) (b0, b1, b2) =
  ( gsub (gmul a1 b2) (gmul a2 b1)
  , gsub (gmul a2 b0) (gmul a0 b2)
  , gsub (gmul a0 b1) (gmul a1 b0)
  )

-- | The CFS rewriting operator: conjugated cross-product.
--
-- This is binary by signature; @compose a b@ ignores any third vertex.
compose :: Cochain -> Cochain -> Cochain
compose a b =
  let (x, y, z) = cross3 a b
  in (gconj x, gconj y, gconj z)

-- | Exact projective ray-equivalence: @a ~ b@ iff @a x b = 0@ (linearly
-- dependent over C). No floating-point tolerance.
projEquiv :: Cochain -> Cochain -> Bool
projEquiv a b =
  let (x, y, z) = cross3 a b
  in gisZero x && gisZero y && gisZero z

isZeroVec :: Cochain -> Bool
isZeroVec (x, y, z) = gisZero x && gisZero y && gisZero z

conjVec :: Cochain -> Cochain
conjVec (x, y, z) = (gconj x, gconj y, gconj z)
