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
    unit $ quietly $ cmd "ghc -o" target "--make" source "-outputdir" outputDir "-hpcdir" outputDir "-j4 -O2 -Wall -Werror -fwarn-unused-do-bind -fno-warn-type-defaults -threaded -rtsopts -fno-ignore-asserts" flags
    unit $ quietly $ cmd "mv -f" (outputDir ++ "/Main.o") (target <.> ".o")
    unit $ quietly $ cmd "mv -f" (outputDir ++ "/Main.hi") (target <.> ".hi")
 
main :: IO ()
main = shakeArgs shakeOptions{shakeFiles="_build"} $ do
    phony "run" $ do
        need ["_build/cfourier"]
        cmd "time _build/cfourier"

    phony "run_hs" $ do
        need ["_build/hsfourier"]
        cmd "time _build/hsfourier +RTS -s"

    phony "view_bench" $ do
        need ["_build/benchmark.html"]
        cmd "firefox _build/benchmark.html"

    "_build/cfourier" %> \out -> do
        need ["cfourier.cc"]
        cmd "g++ -o _build/cfourier cfourier.cc --std=c++11 -O2 -Wall -lOpenCL"

    "_build/hsfourier" %> \out -> buildBin "_build/hsfourier" "hsfourier.hs" "_build/" ""
    "_build/Benchmark" %> \out -> buildBin "_build/Benchmark" "Benchmark.hs" "_build/" ""

    "_build/benchmark.html" %> \out -> do
        need ["_build/Benchmark"]
        cmd "_build/Benchmark --output _build/benchmark.html"
