import Control.Arrow
import System.Random
import Data.List

seed = mkStdGen 0

generateNumbers :: StdGen -> [Int]
generateNumbers = take 1000 . unfoldr (Just . randomR (0, 2^20))

main = putStrLn (show percentage ++ "% of quotients were below 1000")
    where
    (numerators, denominators) = generateNumbers *** generateNumbers $ split seed
    quotients = zipWith div numerators denominators
    percentage = (100 * count (<=1000) quotients `div` length quotients)
    count = (length .) . filter
