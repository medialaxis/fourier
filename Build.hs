import Development.Shake
import Development.Shake.Command
import Development.Shake.FilePath
import Development.Shake.Util

buildBin :: String -- ^ Target.
         -> String -- ^ Main source.
         -> String -- ^ Output directory.
         -> String -- ^ Additional flags.
         -> Action ()
buildBin target source outputDir flags = do
    sources <- getDirectoryFiles "." ["*.hs"]
    need sources
    -- We need to save and restore Main.o and Main.hi when we use ghc with the
    -- -outputdir flag. Otherwise, GHC will use the same output files(Main.o,
    -- Main.hi) for all executables.
    Exit _ <- quietly $ cmd "mv -f" (target <.> ".o") (outputDir ++ "/Main.o")
    Exit _ <- quietly $ cmd "mv -f" (target <.> ".hi") (outputDir ++ "/Main.hi")
    unit $ quietly $ cmd "ghc -o" target "--make" source "-outputdir" outputDir "-hpcdir" outputDir "-j4 -O2 -threaded -rtsopts -fno-ignore-asserts" flags
    unit $ quietly $ cmd "mv -f" (outputDir ++ "/Main.o") (target <.> ".o")
    unit $ quietly $ cmd "mv -f" (outputDir ++ "/Main.hi") (target <.> ".hi")
 
main :: IO ()
main = shakeArgs shakeOptions{shakeFiles="_build"} $ do
    phony "run" $ do
        need ["_build/cfourier"]
        cmd "_build/cfourier"

    "_build/cfourier" %> \out -> do
        need ["cfourier.cc"]
        cmd "g++ -o _build/cfourier cfourier.cc --std=c++11 -O2 -Wall -lOpenCL"
