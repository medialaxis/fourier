import Fourier
import qualified Data.Vector.Unboxed as VU
import Data.Complex

main = do
    let example = VU.replicate (1024*6) ((1 :+ 0) :: Complex Double)
    print $ prop_inverseDft example
