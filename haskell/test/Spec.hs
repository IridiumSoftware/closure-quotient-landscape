{-# LANGUAGE BangPatterns #-}
-- |
-- Cross-implementation validation suite for the CFS rewriting kernel.
--
-- Confirms agreement with @q102_build_exact_v1.py@ on cardinalities for the
-- canonical family (Q24, Q48, Q45, Q51, Q84, Q90, Q102), J fixed-point-free
-- involution, and F-closure terminality on Q48, Q84, Q90, Q102.
module Main (main) where

import qualified Data.Map.Strict as Map
import           System.Exit     (exitFailure, exitSuccess)
import           System.IO       (hFlush, stdout)

import           CFS.Multiway
import           CFS.Topologies

-- Suppress unused-import warning by referring to Map below; this is benign.
_unused :: Map.Map Int Int -> Map.Map Int Int
_unused = id

-- | A single named check producing Bool.
data Check = Check { checkName :: String, checkPass :: Bool, checkDetail :: String }

reportCheck :: Check -> IO ()
reportCheck c = do
  putStrLn $ (if checkPass c then "  PASS  " else "  FAIL  ")
             ++ checkName c
             ++ "  " ++ checkDetail c
  hFlush stdout

main :: IO ()
main = do
  putStrLn "CFS kernel — cross-implementation validation"
  putStrLn "(faithful Haskell port of q102_build_exact_v1.py)"
  putStrLn ""

  putStrLn "--- single-side cardinalities (M(seed)/~) ---"
  let cards =
        [ ("Q_24  (G_0 single side)",        g0Topo,            24, 5)
        , ("Q_42  (adjacent single side)",   adjacentTernary,   42, 5)
        , ("Q_45  (adj+6pads single side)",  q45SeedEdges,      45, 5)
        , ("Q_51  (K_6^3 single side)",      completeTernary 6, 51, 5)
        ]
  cCards <- mapM (\(n, t, e, d) -> singleSideCheck n t e d) cards
  mapM_ reportCheck cCards

  putStrLn ""
  putStrLn "--- C-closed cardinalities ((M ∪ C(M))/~) ---"
  let cclosed =
        [ ("Q_48   (G_0 C-closed)",       g0Topo,            48,  5)
        , ("Q_84   (adjacent C-closed)",  adjacentTernary,   84,  5)
        , ("Q_90   (adj+6pads C-closed)", q45SeedEdges,      90,  5)
        , ("Q_102  (K_6^3 C-closed)",     completeTernary 6, 102, 5)
        ]
  cCC <- mapM (\(n, t, e, d) -> cclosedCheck n t e d) cclosed
  mapM_ reportCheck cCC

  putStrLn ""
  putStrLn "--- J: charge-conjugation involution on Q ---"
  cJ <- mapM jCheck
    [ ("Q_48",  g0Topo)
    , ("Q_84",  adjacentTernary)
    , ("Q_90",  q45SeedEdges)
    , ("Q_102", completeTernary 6)
    ]
  mapM_ reportCheck (concat cJ)

  putStrLn ""
  putStrLn "--- F-closure terminality (F(Q) = Q at depths 4, 5, 6) ---"
  cF <- mapM fCheck
    [ ("Q_48",  g0Topo,            48)
    , ("Q_84",  adjacentTernary,   84)
    , ("Q_90",  q45SeedEdges,      90)
    , ("Q_102", completeTernary 6, 102)
    ]
  mapM_ reportCheck (concat cF)

  putStrLn ""
  let allChecks = cCards ++ cCC ++ concat cJ ++ concat cF
      passed = length (filter checkPass allChecks)
      total  = length allChecks
  putStrLn $ "summary: " ++ show passed ++ " / " ++ show total ++ " checks passed"
  if passed == total then exitSuccess else exitFailure

singleSideCheck :: String -> [Hyperedge] -> Int -> Int -> IO Check
singleSideCheck name seed expected depth = do
  let n = singleSideCardinality seed (canonicalICs 6) depth
  pure $ Check name (n == expected)
           ("got " ++ show n ++ ", expect " ++ show expected)

cclosedCheck :: String -> [Hyperedge] -> Int -> Int -> IO Check
cclosedCheck name seed expected depth = do
  let q  = buildCClosedQuotient seed (canonicalICs 6) depth
      n  = qNcl q
  pure $ Check name (n == expected)
           ("got " ++ show n ++ ", expect " ++ show expected)

jCheck :: (String, [Hyperedge]) -> IO [Check]
jCheck (name, seed) = do
  let q  = buildCClosedQuotient seed (canonicalICs 6) 5
      j  = buildJ q
      iv = jIsInvolution j
      fp = jIsFixedPointFree j
      sz = Map.size j == qNcl q
  pure
    [ Check (name ++ " J total")             sz (show (Map.size j) ++ " / " ++ show (qNcl q))
    , Check (name ++ " J^2 = id")            iv ""
    , Check (name ++ " J fixed-point-free")  fp ""
    ]

fCheck :: (String, [Hyperedge], Int) -> IO [Check]
fCheck (name, seed, expected) = do
  let q = buildCClosedQuotient seed (canonicalICs 6) 5
      cards = [fClosureCardinality q d | d <- [4, 5, 6]]
      ok = all (== expected) cards
  pure
    [ Check (name ++ " F-closure d∈{4,5,6}") ok
        ("|F^d Q| = " ++ show cards ++ " expect " ++ show expected)
    ]

