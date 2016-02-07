module Fourier (
    dft,
    idft,
    prop_inverseDft
    ) where

import qualified Data.Vector.Unboxed as VU
import Data.Complex

type Signal = VU.Vector (Complex Double)

dot :: Complex Double -> Complex Double -> Complex Double
dot a b = a * conjugate b

dotS :: Signal -> Signal -> Complex Double
dotS a b = VU.sum $ VU.zipWith dot a b

minus :: Signal -> Signal -> Signal
minus a b = VU.zipWith (-) a b 

-- | RMS of error signal.
errorS :: Signal -> Signal -> Double
errorS a b = sqrt (realPart (dotS e e)/fromIntegral sz) where
    e = a `minus` b
    sz = VU.length a

i = 0 :+ 1

dft :: Signal -> Signal
dft signal = VU.generate sz go where
    go :: Int -> Complex Double
    go k = sum $ map (f k) [0..sz-1]

    sz :: Int
    sz = VU.length signal

    f :: Int -> Int -> Complex Double
    f k n = (signal VU.! n) * exp (-i*2*pi*fromIntegral k*fromIntegral n/fromIntegral sz)

idft :: Signal -> Signal
idft spectrum = VU.generate sz go where
    go :: Int -> Complex Double
    go n = (sum $ map (f n) [0..sz-1])/fromIntegral sz

    sz :: Int
    sz = VU.length spectrum

    f :: Int -> Int -> Complex Double
    f n k = (spectrum VU.! k) * exp (i*2*pi*fromIntegral k*fromIntegral n/fromIntegral sz)

prop_inverseDft :: Signal -> Double
prop_inverseDft testSignal = errorS testSignal (idft (dft testSignal))
