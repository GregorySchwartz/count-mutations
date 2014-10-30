-- germline_clone_mutation_count_source
-- Mutation Module
-- By G.W. Schwartz

-- Collects all functions pertaining to mutations of codons, including
-- Hamming distances and the like.

module Mutation where

-- Built in
import Data.List
import Data.Maybe
import Control.Applicative

-- Local
import Types
import Translation

-- Checks if a Mutation is four fold redundant
isFourFoldRedundantMutation :: Mutation -> Bool
isFourFoldRedundantMutation x = elem x (mutList list1)
                             || elem x (mutList list2)
                             || elem x (mutList list3)
                             || elem x (mutList list4)
                             || elem x (mutList list5)
                             || elem x (mutList list6)
                             || elem x (mutList list7)
                             || elem x (mutList list8)
  where
    list1 = ["CTT", "CTC", "CTA", "CTG"]
    list2 = ["GTT", "GTC", "GTA", "GTG"]
    list3 = ["TCT", "TCC", "TCA", "TCG"]
    list4 = ["CCT", "CCC", "CCA", "CCG"]
    list5 = ["ACT", "ACC", "ACA", "ACG"]
    list6 = ["GCT", "GCC", "GCA", "GCG"]
    list7 = ["CGT", "CGC", "CGA", "CGG"]
    list8 = ["GGT", "GGC", "GGA", "GGG"]
    mutList a = filter (\(i, j) -> i /= j) $ (,) <$> a <*> a

-- Takes two strings, returns Hamming distance
hamming :: String -> String -> Int
hamming xs ys = length $ filter not $ zipWith (==) xs ys


-- Checks if a pair is actually a mutation
isMutation :: Mutation -> Bool
isMutation (x, y)
    | codon2aa x == codon2aa y = False
    | otherwise                = True

-- Checks if a pair is actually a mutation
isSilentMutation :: Mutation -> Bool
isSilentMutation (x, y)
    | x /= y && codon2aa x == codon2aa y = True
    | otherwise = False

-- Takes a list of mutations and returns the mutations that are
-- actual mutations (the mutation is not just the same character with no gaps.
filterMutStab :: (Mutation -> Bool)
              -> [Mutation]
              -> [Mutation]
filterMutStab isWhat = filter filterRules
  where
    filterRules x    = isWhat x
                    && not (inTuple '-' x)
                    && not (inTuple '.' x)
                    && not (inTuple '~' x)
    inTuple c (x, y) = if c `elem` x || c `elem` y then True else False

-- Return the mutation steps for a mutation
mutation :: Mutation -> Maybe [[(Position, Nucleotide, Bool, Bool)]]
mutation (x, y)
    | hamming x y == 0 = Nothing
    | hamming x y == 2 = Just [ mutSteps (mutPos x y 0 []) x y
                              , mutSteps (reverse . mutPos x y 0 $ []) x y ]
    | hamming x y == 3 = Just [ mutSteps [0, 1, 2] x y
                              , mutSteps [0, 2, 1] x y
                              , mutSteps [1, 0, 2] x y
                              , mutSteps [1, 2, 0] x y
                              , mutSteps [2, 0, 1] x y
                              , mutSteps [2, 1, 0] x y ]
    | otherwise        = Just [mutSteps (mutPos x y 0 []) x y]

-- | A spanning fold that collects the mutation steps from germline to
-- clone codon
-- mutSteps (order of positions mutated) (germline codon) (clone codon)
mutSteps :: [Position] -> Codon -> Codon -> [(Position, Nucleotide, Bool, Bool)]
mutSteps [] _ _     = []
mutSteps (n:ns) x y = ( n
                      , y !! n
                      , isSilentMutation intermediateMutation
                      , isFourFoldRedundantMutation intermediateMutation )
                    : mutSteps ns (changeXToY n x y) y
  where
    intermediateMutation = (x, changeXToY n x y)

-- | Change one nucleotide from xs to ys at position (n - 1) (index at 0)
changeXToY :: Position -> String -> String -> String
changeXToY 0 (_:xs) (y:_) = y:xs
changeXToY n (x:xs) (_:ys) = x : changeXToY (n - 1) xs ys

-- | Determine which positions are mutated in a codon
mutPos :: String -> String -> Position -> [Position] -> [Position]
mutPos _ _ 3 _  = []
mutPos [] _ _ _ = []
mutPos _ [] _ _ = []
mutPos (x:xs) (y:ys) n ns
    | x /= y    = n : mutPos xs ys (n + 1) ns
    | otherwise = mutPos xs ys (n + 1) ns

-- | Find the number of unique synonymous or non-synonymous mutations by
-- nucleotide in each codon while taking into account theoretical
-- intermediate steps from the germline.
uniqueSynonymous :: MutationType
                 -> Bias
                 -> Bool
                 -> CodonMut
                 -> MutCount
                 -> [[Mutation]]
                 -> Int
uniqueSynonymous mutType bias fourBool codonMut mutCount = sum
                                                         . map ( length
                                                         . getMutationCount )
  where
    getMutationCount   = nub -- Unique mutations
                       . mutCountFrequent
                       . filter (\(_, _, sil, four) -> biasValue mutType sil
                                                    && isFourFold fourBool four)
                       . concatMap (mutBias . fromJust)
                       . filter isJust
                       . map (mutatedCodon codonMut)  -- Only codons with n muts
                       . map mutation
    mutCountFrequent = concat
                     . filter (\x -> length x >= mutCount)
                     . group
                     . sort
    mutBias []       = []
    mutBias xs       = xs !! (fromJust $ biasResult xs)
    biasResult xs    = elemIndex (biasIndex mutType bias xs) (map sumMut xs)
    biasIndex Silent Silent           = maximum . map sumMut
    biasIndex Silent Replacement      = minimum . map sumMut
    biasIndex Replacement Silent      = minimum . map sumMut
    biasIndex Replacement Replacement = maximum . map sumMut
    sumMut          = sum
                    . map ( \(_, _, sil, _)
                         -> if ((biasValue mutType) sil) then 1 else 0 )
    biasValue Silent       = id
    biasValue Replacement  = not
    isFourFold True x      = x  -- Check four fold redundancy
    isFourFold False _     = True  -- Include it all if we don't care about it
    mutatedCodon 0 xs      = xs
    mutatedCodon _ Nothing = Nothing
    mutatedCodon 1 (Just xs)
        | length xs == 1   = Just xs
        | otherwise        = Nothing
    mutatedCodon 2 (Just xs)
        | length xs == 2   = Just xs
        | otherwise        = Nothing
    mutatedCodon 3 (Just xs)
        | length xs == 6   = Just xs
        | otherwise        = Nothing

-- Return the mutations if they exist between the germline and a clone
-- (unused in this algorithm, only here for completionist reasons).
countMutations :: Germline
               -> Clone
               -> [(Position, Mutation)]
countMutations germline clone = mut germline clone
  where
    mut x y         = zip [1..] . zip x $ y
