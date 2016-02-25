import Fourier
import Criterion.Main
import qualified Data.Vector.Unboxed as VU
import Data.Complex

main :: IO ()
main = defaultMain
    [ bgroup "hsfourier"
        [ bench "prop_inverseDft" $
            let example = VU.replicate (1024*6) ((1 :+ 0) :: Complex Double)
            in whnf prop_inverseDft example
        ]
    ]
