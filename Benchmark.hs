import Fourier
import Criterion.Main
import qualified Data.Vector.Unboxed as VU
import Data.Complex

fib :: Integer -> Integer
fib m | m < 0     = error "negative!"
      | otherwise = go m
  where
    go 0 = 0
    go 1 = 1
    go n = go (n-1) + go (n-2)

main :: IO ()
main = defaultMain
    [ bgroup "fib"
        [ bench "1" $ whnf fib 1
        , bench "9" $ whnf fib 9
        ]
    , bgroup "hsfourier"
        [ bench "prop_inverseDft" $
            let example = VU.replicate (1024*6) ((1 :+ 0) :: Complex Double)
            in whnf prop_inverseDft example
        ]
    ]
