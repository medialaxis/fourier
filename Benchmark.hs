-- import Fourier
import Criterion.Main

fib :: Integer -> Integer
fib m | m < 0     = error "negative!"
      | otherwise = go m
  where
    go 0 = 0
    go 1 = 1
    go n = go (n-1) + go (n-2)

--    let example = VU.replicate (1024*6) ((1 :+ 0) :: Complex Double)
--    print $ prop_inverseDft example

main :: IO ()
main = defaultMain
    [ bgroup "hsfourier"
        [ bench "1" $ whnf fib 1
        , bench "9" $ whnf fib 9
        ]
    ]
