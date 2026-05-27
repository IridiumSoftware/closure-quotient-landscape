{-# LANGUAGE BangPatterns #-}
-- |
-- Module      : CFS.Multiway
-- Description : Ternary-hyperedge multiway rewriting, gauge quotient, J doubling.
--
-- This module makes the ternary structure of the CFS rewriting graph
-- type-explicit:
--
-- > type Hyperedge = (Vertex, Vertex, Vertex)
--
-- Each step takes a hyperedge @(v1, v2, v3)@, applies the /binary/ operator
-- 'CFS.Cochain.compose' to @(psi v1, psi v2)@ to produce a new cochain @w@,
-- introduces a new vertex carrying @w@, and schedules /three/ new hyperedges
-- propagating @w@, @v1@, @v2@, and @v3@. The third slot is /not/ an input
-- to the operator but is genuinely part of the rewriting graph's structure:
-- it shapes which hyperedges fire at depth d+1 and thus participates in the
-- gauge-quotient cardinality |Q|.
--
-- Faithful Haskell port of @q102_build_exact_v1.py@.
module CFS.Multiway
  ( Vertex
  , Hyperedge
  , Depth
  , buildMultiway
  , Quotient(..)
  , Tier(..)
  , ClOrigin(..)
  , buildCClosedQuotient
  , singleSideCardinality
    -- * J: charge-conjugation involution on Q
  , JMap
  , buildJ
  , jIsInvolution
  , jIsFixedPointFree
    -- * F-closure terminality
  , fClosureCardinality
  ) where

import           CFS.Cochain
import qualified Data.Map.Strict as Map
import           Data.Map.Strict (Map)
import qualified Data.Set        as Set
import           Data.Set        (Set)

type Vertex = Int

-- | A ternary hyperedge: three vertex indices.
type Hyperedge = (Vertex, Vertex, Vertex)

type Depth = Int

-- | Faithful port of @build_multiway@ from q102_build_exact_v1.py.
--
-- Each step takes the hyperedges scheduled at the current depth, applies the
-- /binary/ 'compose' to the first two vertices' cochains, introduces a new
-- vertex carrying the result, and schedules three new hyperedges for the
-- next depth. Composition is /cached/ on the @(v1, v2)@ pair, so two
-- hyperedges sharing the same first two slots reuse the same new vertex.
buildMultiway
  :: [Hyperedge]                     -- ^ seed hyperedges
  -> Map Vertex Cochain               -- ^ initial cochains psi(v) for seed vertices
  -> Depth                            -- ^ rewriting depth
  -> (Map Vertex Cochain, [(Depth, Vertex, Vertex, Vertex)])
buildMultiway seedEdges psi0 depth =
  let startNext = case Map.lookupMax psi0 of
        Just (k, _) -> k + 1
        Nothing     -> 0
      seeded    = [(0 :: Depth, v1, v2, v3) | (v1, v2, v3) <- seedEdges]
      (psiF, _, edgesF) = sliceLoop startNext psi0 Map.empty seeded 0
  in (psiF, edgesF)
  where
    sliceLoop
      :: Vertex
      -> Map Vertex Cochain
      -> Map (Vertex, Vertex) Vertex
      -> [(Depth, Vertex, Vertex, Vertex)]
      -> Depth
      -> ( Map Vertex Cochain
         , Map (Vertex, Vertex) Vertex
         , [(Depth, Vertex, Vertex, Vertex)]
         )
    sliceLoop !nxt !psi !cache !edges !d
      | d >= depth = (psi, cache, edges)
      | otherwise  =
          let here = [e | e@(de, _, _, _) <- edges, de == d]
              (nxt', psi', cache', newE) =
                foldl stepEdge (nxt, psi, cache, []) here
          in sliceLoop nxt' psi' cache' (edges ++ reverse newE) (d + 1)

    stepEdge
      :: ( Vertex
         , Map Vertex Cochain
         , Map (Vertex, Vertex) Vertex
         , [(Depth, Vertex, Vertex, Vertex)]
         )
      -> (Depth, Vertex, Vertex, Vertex)
      -> ( Vertex
         , Map Vertex Cochain
         , Map (Vertex, Vertex) Vertex
         , [(Depth, Vertex, Vertex, Vertex)]
         )
    stepEdge (nxt, psi, cache, acc) (d, v1, v2, v3) =
      let key = (v1, v2)
      in case Map.lookup key cache of
           Just w ->
             let acc' = (d + 1, w, v1, v2)
                      : (d + 1, w, v1, v3)
                      : (d + 1, w, v2, v3)
                      : acc
             in (nxt, psi, cache, acc')
           Nothing ->
             case (Map.lookup v1 psi, Map.lookup v2 psi) of
               (Just a, Just b) ->
                 let w = compose a b
                 in if isZeroVec w
                      then (nxt, psi, cache, acc)
                      else
                        let psi'   = Map.insert nxt w psi
                            cache' = Map.insert key nxt cache
                            acc'   = (d + 1, nxt, v1, v2)
                                   : (d + 1, nxt, v1, v3)
                                   : (d + 1, nxt, v2, v3)
                                   : acc
                        in (nxt + 1, psi', cache', acc')
               _ -> (nxt, psi, cache, acc)

-- | Tier of a quotient cluster: A (origin IC), B (depth-1 child), C (deeper).
data Tier = TierA | TierB | TierC
  deriving (Eq, Ord, Show)

data ClOrigin = OrigOnly | ConjOnly | Both
  deriving (Eq, Ord, Show)

-- | The exact C-closed quotient Q = (M(seed) ∪ C(M(seed)))/~.
data Quotient = Quotient
  { qNcl      :: !Int                 -- ^ cluster count |Q|
  , qPsi      :: !(Map Int Cochain)   -- ^ representative cochain per cluster
  , qHe       :: !(Set Hyperedge)     -- ^ induced hyperedges on clusters
  , qTier     :: !(Map Int Tier)
  , qClOrigin :: !(Map Int ClOrigin)
  } deriving (Show)

-- | Faithful port of @build_c_closed_quotient@: build orig + J-conjugated
-- multiway, merge, exact-ray-cluster, return Q.
buildCClosedQuotient
  :: [Hyperedge]              -- ^ seed topology
  -> Map Vertex Cochain        -- ^ canonical Gaussian-integer ICs
  -> Depth                     -- ^ rewriting depth
  -> Quotient
buildCClosedQuotient seedEdges ics depth =
  let (psiO, edgesO)    = buildMultiway seedEdges ics depth
      icsConj           = Map.map conjVec ics
      (psiC, edgesCRaw) = buildMultiway seedEdges icsConj depth
      offset            = case Map.lookupMax psiO of
        Just (k, _) -> k + 1
        Nothing     -> 0
      psiAll  = Map.union psiO (Map.mapKeysMonotonic (+ offset) psiC)
      origIds = Map.keysSet psiO
      edgesC  = [(d, a + offset, b + offset, c + offset)
                | (d, a, b, c) <- edgesCRaw]

      verts = [(v, p) | (v, p) <- Map.toAscList psiAll, not (isZeroVec p)]
      (vidToCid, reps) = clusterByRay verts
      nCl = length reps

      icKeys = Map.keysSet ics
      origCids = Set.fromList
        [ c | v <- Set.toList icKeys, Just c <- [Map.lookup v vidToCid] ]

      -- Depth-1 vertices from BOTH the orig and conj sides (conj side
      -- vertices are offset).
      depth1Vids =
        Set.fromList [ w | (1, w, _, _) <- edgesO ]
          `Set.union`
        Set.fromList [ w + offset | (1, w, _, _) <- edgesCRaw ]
      gen1 = Set.fromList
        [ c | v <- Set.toList depth1Vids
            , Just c <- [Map.lookup v vidToCid]
            , not (Set.member c origCids)
        ]

      tierOf c
        | Set.member c origCids = TierA
        | Set.member c gen1     = TierB
        | otherwise             = TierC

      -- Origin breakdown per cluster: which sides contributed members
      bumpOrigin :: Vertex -> Map Int (Bool, Bool) -> Int -> Map Int (Bool, Bool)
      bumpOrigin v m c =
        let isOrig = Set.member v origIds
            (o, j) = Map.findWithDefault (False, False) c m
        in if isOrig
             then Map.insert c (True, j) m
             else Map.insert c (o, True) m
      sourcePerCl =
        foldl (\m (v, _) -> case Map.lookup v vidToCid of
                              Nothing -> m
                              Just c  -> bumpOrigin v m c)
              Map.empty verts
      clOriginOf c = case Map.findWithDefault (False, False) c sourcePerCl of
        (True, True)  -> Both
        (True, False) -> OrigOnly
        _             -> ConjOnly

      he = Set.fromList
        [ (c1, c2, c3)
        | (_, v1, v2, v3) <- edgesO ++ edgesC
        , Just c1 <- [Map.lookup v1 vidToCid]
        , Just c2 <- [Map.lookup v2 vidToCid]
        , Just c3 <- [Map.lookup v3 vidToCid]
        ]
  in Quotient
       { qNcl      = nCl
       , qPsi      = Map.fromList (zip [0 ..] reps)
       , qHe       = he
       , qTier     = Map.fromList [(c, tierOf c) | c <- [0 .. nCl - 1]]
       , qClOrigin = Map.fromList [(c, clOriginOf c) | c <- [0 .. nCl - 1]]
       }

-- | Single-side cardinality: |M(seed)/~| with no C-closure.
singleSideCardinality
  :: [Hyperedge]
  -> Map Vertex Cochain
  -> Depth
  -> Int
singleSideCardinality seedEdges ics depth =
  let (psi, _) = buildMultiway seedEdges ics depth
      verts    = [(v, p) | (v, p) <- Map.toAscList psi, not (isZeroVec p)]
      (_, reps) = clusterByRay verts
  in length reps

-- | Cluster cochains by exact projective ray equivalence. Returns the
-- vertex-to-cluster map and the ordered list of representatives.
clusterByRay :: [(Vertex, Cochain)] -> (Map Vertex Int, [Cochain])
clusterByRay = go Map.empty [] 0
  where
    go !v2c !reps !_n [] = (v2c, reverse reps)
    go !v2c !reps !n ((v, p) : rest) =
      case findRay reps p of
        Just c  -> go (Map.insert v c v2c) reps n rest
        Nothing -> go (Map.insert v n v2c) (p : reps) (n + 1) rest

    -- 'reps' is in reverse order (newest first); we search in original order
    -- so existing clusters bind to the lowest index.
    findRay reps p =
      let ordered = zip [0 :: Int ..] (reverse reps)
      in case [c | (c, r) <- ordered, projEquiv p r] of
           (c : _) -> Just c
           []      -> Nothing

-- | Charge-conjugation involution J on the cluster set.
type JMap = Map Int Int

buildJ :: Quotient -> JMap
buildJ q =
  let reps = qPsi q
      ns   = [0 .. qNcl q - 1]
      go c = do
        r <- Map.lookup c reps
        let cr = conjVec r
        case [c' | c' <- ns
                , Just r' <- [Map.lookup c' reps]
                , projEquiv cr r'] of
          (c' : _) -> Just (c, c')
          []       -> Nothing
  in Map.fromList [ (c, c') | c <- ns, Just (_, c') <- [go c] ]

-- | J should square to the identity on the cluster set.
jIsInvolution :: JMap -> Bool
jIsInvolution j = all step (Map.toList j)
  where
    step (c, c') = case Map.lookup c' j of
      Just c'' -> c'' == c
      Nothing  -> False

-- | J should have no fixed points (charge-conjugation flips every cluster).
jIsFixedPointFree :: JMap -> Bool
jIsFixedPointFree = all (\(c, c') -> c /= c') . Map.toList

-- | Re-seed the multiway with Q's cluster cochains as initial cochains and
-- Q's induced hyperedges as the topology, then count clusters at the given
-- depth. Equality with @qNcl q@ at depths 4, 5, 6 means Q is F-closed (i.e.,
-- doubly terminal alongside C-closure).
fClosureCardinality :: Quotient -> Depth -> Int
fClosureCardinality q d =
  let psi0  = qPsi q
      edges = Set.toAscList (qHe q)
      (psi', _) = buildMultiway edges psi0 d
      verts     = [(v, p) | (v, p) <- Map.toAscList psi', not (isZeroVec p)]
      (_, reps) = clusterByRay verts
  in length reps
