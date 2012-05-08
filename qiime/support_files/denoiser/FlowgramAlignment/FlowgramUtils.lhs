> module FlowgramUtils where

> import Data.Array
> import System.IO.Unsafe (unsafePerformIO)

> type FlowSignalDistrib = Array (Int,Int) Float

> type Flowgram  = [Float]
> type Flow = (Float, Float, Float, Float)
> type Flowgram_tup = [Flow]
> type FlowOrder = [Char]


> parseSFFfile:: [String] ->[Flowgram]
> parseSFFfile [] = []
> parseSFFfile (x:xs)
>     | x == ""   = parseSFFfile xs
>     | otherwise =  case head (words x) of
>                       "Flowgram:" -> map read (tail (words x)) : parseSFFfile xs 
>                       otherwise -> parseSFFfile xs

> -- beware of side effect!
> flow_sig_dist:: String -> FlowSignalDistrib
> flow_sig_dist filename = unsafePerformIO $ readFlowSignalDistrib 1000 filename

> readFlowSignalDistrib:: Int -> String -> IO FlowSignalDistrib
> readFlowSignalDistrib bins file = do
>                                     s <- readFile file
>                                     let xs = lines s
>                                         n  = length xs  `div` (bins+1) 
>                                         arr= listArray ((0,0),((n-1),(bins-1)))
>                                                        (splitSignalDistrib bins xs)
>                                     return arr

> splitSignalDistrib:: Int -> [String] -> [Float]
> splitSignalDistrib bins [] = []
> splitSignalDistrib bins (x:xs) = ((map read h)::[Float]) ++ splitSignalDistrib bins t
>                    where (h,t) = (take bins xs, drop bins xs)
>                          -- Chris' .dat files contain one extra value (x above) for each signal  
>                          -- Has been used in an old version. Can be ignored


> tupelize_flow:: Flowgram-> Flowgram_tup
> tupelize_flow [] = []
> tupelize_flow (a:b:c:d:xs) = (a,b,c,d): tupelize_flow xs
> tupelize_flow (x:xs) = error ("Flowgram length is not a multiple of 4!")

> def_fo:: FlowOrder
> def_fo = "TCAG"

> flow_to_seq:: FlowOrder -> Flowgram -> [Char]
> flow_to_seq fo flow = concat [ take n (repeat nuc)|  (nuc, signal) <- zip (concat (repeat fo)) flow, let n = round' signal]

> round':: Float -> Int
> round' f = floor (f + 0.5)

