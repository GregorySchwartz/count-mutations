-- count-mutations
-- Main
-- By G.W. Schwartz

-- Takes a fasta file in the format of ">>Germline header\nGermline
-- sequence\n>Mutant header\nMutant sequence...etc" and generates the
-- mutation counts of the clones from each listed germline in the file per
-- position in a dataframe type format. Contains flags to bias for or
-- against silent or replacement mutations, along with only including
-- codons with certain number of mutations.

-- Compiling instructions:
-- Requires the Haskell Platform (sudo apt-get install haskell-platform)
-- Requires Options.Applicative (cabal update && cabal install
-- optparse-applicative)
-- Requires Data.List.Split (cabal update && cabal install split)
-- Compiling command: ghc --make /path/to/Main.hs -o ./mutation_count -O2

-- Built in
import qualified Data.Map as M

-- Cabal
import Options.Applicative

-- Local
import Types
import FastaRead
import MutationCount

-- | Command line arguments
data Options = Options { inputMutType  :: MutationType
                       , inputBias     :: Bias
                       , inputCodonMut :: CodonMut
                       , inputMutCount :: MutCount
                       , inputLabel    :: String
                       , removeN       :: Bool
                       , onlyFourFold  :: Bool
                       , input         :: String
                       , inputFasta    :: String
                       , output        :: String
                       }

-- | Command line options
options :: Parser Options
options = Options
      <$> option auto
          ( long "input-mut-type"
         <> short 't'
         <> metavar "[Silent]|Replacement"
         <> value Silent
         <> help "The type of mutation to be counted" )
      <*> option auto
          ( long "input-bias"
         <> short 'b'
         <> metavar "[Silent]|Replacement"
         <> value Silent
         <> help "The type of mutation to bias for when calculating\
                 \ intermediate codons" )
      <*> option auto
          ( long "input-codon-mut"
         <> short 'c'
         <> metavar "[0]|1|2|3"
         <> value 0
         <> help "Only count mutations from codons with this many mutations\
                 \ (0 is the same as include all codons)" )
      <*> option auto
          ( long "input-mut-count"
         <> short 'm'
         <> metavar "[1]|2|3|..."
         <> value 1
         <> help "Only count a unique mutation if it appears this many\
                 \ or more times" )
      <*> strOption
          ( long "input-label"
         <> short 'l'
         <> metavar "[count]|STRING"
         <> value "count"
         <> help "The label for the data (usually the dataset)" )
      <*> switch
          ( long "remove-N"
         <> short 'N'
         <> help "Whether to remove N or n in the sequence" )
      <*> switch
          ( long "four-fold-redundant"
         <> short 'f'
         <> help "Whether to only count mutations if they are\
                 \ four fold redundant mutations" )
      <*> strOption
          ( long "input"
         <> short 'I'
         <> metavar "STRING STRING"
         <> value ""
         <> help "Two sequences separated by spaces in order to find\
                 \ the mutations between them" )
      <*> strOption
          ( long "input-fasta"
         <> short 'i'
         <> metavar "FILE"
         <> value ""
         <> help "The name of the input fasta file" )
      <*> strOption
          ( long "output"
         <> short 'o'
         <> metavar "FILE"
         <> value ""
         <> help "The name of the output file with the counts" )

-- | Get the input file in the correct format
getInput :: Options -> IO String
getInput opts = do
    if (not . null . inputFasta $ opts)
        then readFile . inputFasta $ opts
        else
            if (not . null . input $ opts)
                then return (">>SEQ1\n" ++ head seqs ++ "\n>SEQ2\n" ++ last seqs)
                else getContents
  where
    seqs = words . input $ opts

mutationCounts :: Options -> IO ()
mutationCounts opts = do
    contents     <- getInput opts
    let mutType  = inputMutType opts
        bias     = inputBias opts
        codonMut = inputCodonMut opts
        mutCount = inputMutCount opts
        label    = inputLabel opts

    -- Get rid of carriages
        contentsNoCarriages  = filter (/= '\r') $ contents
    -- No newlines in sequence
        contentsNoNewlines  = joinSeq (removeN opts) contentsNoCarriages

        cloneMap  = generateCloneMap contentsNoNewlines

    -- Generate the germline to clone map
        cloneMutMap         = generateCloneMutMap cloneMap

    -- Generate the germline to clone map separated by clones
        positionCloneMap   = generatePositionCloneMap cloneMutMap

    -- Generate the position to clone map
        combinedCloneMutMap = M.unionsWith (++)
                            . map snd
                            . M.toAscList
                            $ cloneMutMap
    -- Get the string to save
        outputData          = printMutCounts label
                                            mutType
                                            bias
                                            (onlyFourFold opts)
                                            codonMut
                                            mutCount
                                            combinedCloneMutMap
                                            positionCloneMap

    -- Save the counts
    if (null . output $ opts)
        then putStrLn outputData
        else writeFile (output opts) outputData

main :: IO ()
main = execParser opts >>= mutationCounts
  where
    opts = info (helper <*> options)
      ( fullDesc
     <> progDesc "Return the mutation counts with certain biases from the \
                 \ germline to the mutants within clones"
     <> header "count-mutations, Gregory W. Schwartz" )
