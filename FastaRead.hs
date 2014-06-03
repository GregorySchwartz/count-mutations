-- germline_clone_mutation_count_source
-- FastaRead Module
-- By G.W. Schwartz

-- Collects all functions pertaining to the reading of an input fasta in
-- the format of ">>Germline header\nGermline sequence\n>Mutant
-- header\nMutant sequence...etc".

module FastaRead where

-- Built in
import Data.List
import qualified Data.Map as M

-- Cabal
import qualified Data.List.Split as Split

-- Local
import Types

-- Takes a fasta file string and removes newlines in the sequences to make
-- this compatible with the fasta parser. The lineCompress function should
-- get rid of any extra newlines messing with the code.
joinSeq :: Bool -> String -> String
joinSeq removeNFlag = lineCompress
                    . tail
                    . concat
                    . map newEntry
                    . filter (/= "")
                    . Split.splitOn ">>"
  where
    newEntry x             = if elem '>' x then cloneEntry x else germEntry x
    germEntry x            = newGerm x
    cloneEntry x           = newGerm (germline x)
                          ++ concat (map newClone . filter (/= "") . clone $ x)
    newGerm x
        | seq x /= ""      = "\n>>" ++ (header x) ++ "\n" ++ (seq x)
        | otherwise        = ""
    newClone x
        | seq x /= ""      = "\n>" ++ (header x) ++ "\n" ++ (seq x)
        | otherwise        = ""
    germline               = head . Split.splitOn ">"
    clone                  = tail . Split.splitOn ">"
    header                 = head . lines
    seq                    = map (removeN removeNFlag) . concat . tail . lines
    lineCompress []        = []
    lineCompress ('\n':xs) = '\n' : (lineCompress $ dropWhile (== '\n') xs)
    lineCompress (x:xs)    = x : (lineCompress xs)
    removeN False x        = x
    removeN True 'N'       = '-'
    removeN True 'n'       = '-'
    removeN True x         = x

-- Takes a fasta file string of the format ">>[Germline header]\n[Germline
-- sequence]\n>[Mutant header]\n[Mutant
-- sequence]\n>[MH]\n[MS]...\n>>[GH]\n[GS]\n>[MH]\n[MS]...etc" and returns
-- a CloneMap in order to generate the basic building block for the
-- mutation counting.  Note: Several repeating germlines, so they need
-- a unique identifier (an integer in this case). Also removes empty clones
-- (only a germline appears).
generateCloneMap :: String -> CloneMap
generateCloneMap = M.fromList . map getCodonSplit . getSequences
  where
    getCodonSplit ((x1, x2), xs)   = ((x1, fullCodon . Split.chunksOf 3 $ x2),
                                      map (fullCodon . Split.chunksOf 3) xs)
    fullCodon                      = filter ((==3) . length)
    getSequences                   = map filterAssocList . assocList
    filterAssocList ((x1, x2), xs) = ((x1, last . filterHeaders $ x2),
                                      filterHeaders xs)
    filterHeaders                  = filter (\x -> head x /= '>')
    assocList                      = map assocMap . germlineSplit
    assocMap (x, y)                = ((x, germline y), clones y)
    germlineSplit                  = zip [0..]
                                   . filter (\x -> elem '>' x)
                                   . Split.splitOn ">>"
    germline                       = take 2 . lines
    clones                         = drop 2 . lines
