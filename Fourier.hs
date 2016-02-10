module Fourier (
    dft,
    idft,
    prop_inverseDft
    ) where

import qualified Data.Vector.Unboxed as VU
import Data.Complex

scale :: Num a => a -> Complex a -> Complex a
scale x (a :+ b) = (x*a :+ x*b)

type Signal = VU.Vector (Complex Double)

dot :: Complex Double -> Complex Double -> Complex Double
dot a b = a * conjugate b

dotS :: Signal -> Signal -> Complex Double
dotS a b = VU.sum $ VU.zipWith dot a b

minus :: Signal -> Signal -> Signal
minus a b = VU.zipWith (-) a b 

-- | RMS of error signal.
errorS :: Signal -> Signal -> Double
errorS a b = sqrt ((realPart (dotS e e))/fromIntegral sz) where
    e = a `minus` b
    sz = VU.length a

i = 0 :+ 1

dft :: Signal -> Signal
dft signal = VU.generate sz $ \k -> VU.sum $ VU.imap (f k) signal where
    f :: Int -> Int -> Complex Double -> Complex Double
    f k n sn = sn * exp (scale (-2*pi*fromIntegral k*fromIntegral n/fromIntegral sz) i)

    sz :: Int
    sz = VU.length signal

idft :: Signal -> Signal
idft spectrum = VU.generate sz $ \n -> scale (1/fromIntegral sz) (VU.sum $ VU.imap (f n) spectrum) where
    f :: Int -> Int -> Complex Double -> Complex Double
    f n k sn = sn * exp (scale (2*pi*fromIntegral k*fromIntegral n/fromIntegral sz) i)

    sz :: Int
    sz = VU.length spectrum

prop_inverseDft :: Signal -> Double
prop_inverseDft testSignal = errorS testSignal (idft (dft testSignal))
