-- germline_clone_mutation_count_source
-- Types Module
-- By G.W. Schwartz

-- Collects all application specific types

module Types where

-- Built in
import qualified Data.Map as M

-- Algebraic
data MutationType = Silent | Replacement deriving (Eq, Show, Read)

-- Basic
type ID         = Int
type Nucleotide = Char
type Codon      = String
type Clone      = [Codon]
type Germline   = [Codon]
type Position   = Int
type CodonMut   = Int
type Bias       = MutationType
type MutCount   = Int

-- Advanced
type Mutation         = (Codon, Codon)
type CloneMap         = M.Map (ID, Germline) [Clone]
type MutationMap      = M.Map Position [Mutation]
type CloneMutMap      = M.Map (ID, Germline) MutationMap
-- | Same as the Mutation map but separates the values by clones
type PositionCloneMap = M.Map Position [[Mutation]]
