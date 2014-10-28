-- MutationCount Module
-- By G.W. Schwartz

-- Collects all functions pertaining to the counting of mutations called
-- from Main

module MutationCount where

-- Built in
import Data.List
import Data.Char
import Data.Maybe
import qualified Data.Map as M

-- Local
import Types
import Mutation

-- Join together mutation lists into a MutationMap
joinMutations :: [[(Position, Mutation)]] -> MutationMap
joinMutations = groupedMutations
  where
    groupedMutations = M.fromListWith (++) . map (\(x, y) -> (x, [y])) . concat

-- Generate a CloneMutMap which will then be printed to save files
generateCloneMutMap :: CloneMap -> CloneMutMap
generateCloneMutMap = M.mapWithKey gatherMutations
  where
    gatherMutations k xs = joinMutations . map (countMutations (snd k)) $ xs

-- Generate the position to clone map that is separated by clones
generatePositionCloneMap :: CloneMutMap -> PositionCloneMap
generatePositionCloneMap = joinSilentMutations
                         . M.map (M.map (wrapClone . realCodon))
  where
    joinSilentMutations = M.unionsWith (++) . map snd . M.toAscList
    wrapClone xs        = [xs]
    realCodon           = filterMutStab (\_ -> True)

-- Return the results of the mutation or stable counts as a string
printMutCounts :: String
               -> MutationType
               -> Bias
               -> Bool
               -> CodonMut
               -> MutCount
               -> MutationMap
               -> PositionCloneMap
               -> String
printMutCounts label
               mutType
               bias
               fourBool
               codonMut
               mutCount
               mutationMap
               mutTypeMap = header ++ body
  where
    header           = "label,position,count,weight\n"
    body             = unlines
                     . map mapLine
                     . M.toAscList
                     $ mutationMap
    mapLine (x, xs) = label
                   ++ ","
                   ++ show x
                   ++ ","
                   ++ (show . trueMutCount x $ mutTypeMap)
                   ++ ","
                   ++ (show . length . filterMutStab (\_ -> True) $ xs)
    trueMutCount p  = uniqueSynonymous mutType bias fourBool codonMut mutCount
                    . fromMaybe [[]]
                    . M.lookup p
