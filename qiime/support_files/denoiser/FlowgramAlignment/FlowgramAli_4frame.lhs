> module Main where

> import Data.Array
> import Data.List
> import FlowgramUtils
> import ADPCombinators
> import System.Environment(getArgs)
> import System.IO
> import System.Exit(exitWith, ExitCode(..))
> import Text.Printf(printf)

The signature:

> data FG_Alignment =
>                  Nil                                                             |
>                  D  Float Float Float Float FG_Alignment                         |
>                  I                          FG_Alignment Float Float Float Float |
>                  R  Float                   FG_Alignment                   Float |
>                  Td Subword  ()                                                  |
>                  Ti  () Subword            
>                                            deriving (Eq, Show)

Algebra type:

> type FG_Algebra alphabet answer = (
>   answer,                                                            -- nil
>   alphabet -> alphabet -> alphabet -> alphabet -> answer -> answer,  -- d      deletion
>   answer -> alphabet -> alphabet -> alphabet -> alphabet -> answer,  -- i      insertion
>   alphabet -> answer -> alphabet                         -> answer,  -- r      replacement/match
>   Subword  -> ()                                         -> answer,  -- td     terminal deletion
>               ()     -> Subword                          -> answer,  -- ti     terminal insertion
>   [answer]                                               -> [answer] -- h      choice function
>   )

Enumeration algebra:

> enum :: FG_Algebra Float FG_Alignment
> enum = (nil, d, i, r, td, ti, h) where
>    nil = Nil
>    d   = D
>    i   = I
>    r   = R
>    td  = Td
>    ti  = Ti
>    h   = id

Pretty printing algebra:

> prettyprint ::  Array Int Float -> FG_Algebra Float (String, String)
> prettyprint inp = (nil, d, i, r, td, ti, h) where
>   nil          = ("","")
>   d  a b c d (l,r)         = (show a ++ " "++show b ++" "++ show c ++ " "++show d ++ " " ++ l, gap ++ r)
>   i          (l,r) a b c d = (gap ++ l, show d ++" "++ show c ++" "++ show b ++" "++ show a ++ " " ++ r)
>   r  x       (l,r)       y = (show x ++ " " ++ l, show y ++ " " ++ r)
>   td (i,j) ()              = (show[  (inp!k) | k <- [i+1..j]], "")
>   ti () (i,j)              = ("", concat [ " "++ show(inp!k) | k <- [i+1..j]])
>   h            = id
>   gap          = " -  -  -  - "

Pretty printing algebra:
(Debugging algebra)

> prettyprint' ::  FlowSignalDistrib -> Array Int Float -> FG_Algebra Float (String, String, [Float])
> prettyprint' arr inp = (nil, d, i, r, td, ti, h) where
>   nil          = ("","", [])
>   d  a b c d (l,r,s)         = (show a ++ " "++show b ++" "++ show c ++ " "++show d ++ " " ++ l, gap ++ r,
>                                 [15.0, 15, 15, 15] ++ s)
>   i          (l,r,s) a b c d = (gap ++ l, show d ++" "++ show c ++" "++ show b ++" "++ show a ++ " " ++ r,
>                                 [15.0, 15, 15, 15] ++ s)
>   r  x       (l,r,s)       y = (show x ++ " " ++ l, show y ++ " " ++ r, scr:s)
>       where scr = case round' x > 9 of 
>                   True  ->  arr!(       9, floor( (min y 9.99) * 100))
>                   False ->  arr!(round' x, floor( (min y 9.99) * 100))
>   td (i,j) ()              = (show[  (inp!k) | k <- [i+1..j]], "",[])
>   ti () (i,j)              = ("", concat [ " "++show(inp!k) | k <- [i+1..j]],[])
>   h            = id
>   gap          = " -  -  -  - "


Counting algebra:

> count :: FG_Algebra Float Int
> count = (nil, d, i, r, td, ti, h) where
>    nil     = 1
>    d _ _ _ _ s   = s
>    i s _ _ _ _   = s
>    r a s b = s
>    td  r ()  = 1
>    ti () r   = 1
>    h []    = []
>    h l     = [sum l]

score algebra:

> score :: FlowSignalDistrib -> FG_Algebra Float Float
> score arr = (nil, d, i, r, td, ti, h) where
>    nil     = 0
>    d _ _ _ _ s         = s + gap * 4
>    i         s _ _ _ _ = s + gap * 4
>    r a       s       b = case round' a > 9 of 
>                          True  -> s +  arr!(       9, floor( (min b 9.99) * 100))
>                          False -> s +  arr!(round' a, floor( (min b 9.99) * 100))
>    td  r ()  = 0
>    ti () r   = 0
>    h []      = []
>    h l       = [minimum l]

>    gap = 15

mismatch algebra:
Minimizes mismatches between translated nucleotide sequences.

> mismatch ::  FG_Algebra Float Int
> mismatch = (nil, d, i, r, td, ti, h) where
>    nil     = 0
>    d a b c d s         = s + round' a + round' b + round' c + round' d
>    i         s a b c d = s + round' a + round' b + round' c + round' d
>    r a       s       b = s + abs(round' a - round' b)
>                                
>    td  r ()  = 0   -- only core alignment counts
>    ti () r   = 0
>    h []      = []
>    h l       = [minimum l]

(mismatch x aligned seq length) algebra:
Second element of tuple counts the core alignments length of 
       the translated nucleotide seqs

> mismatch_seqlen ::  FG_Algebra Float (Int,Int)
> mismatch_seqlen = (nil, d, i, r, td, ti, h) where
>    nil     = (0,0)
>    d a b c d (s,t)         = (s+x, t+x) 
>        where x = round' a + round' b + round' c + round' d
>    i         (s,t) a b c d = (s+x, t+x)
>        where x = round' a + round' b + round' c + round' d
>    r a       (s,t)       b = (s + abs(round' a - round' b), t+max (round' a) (round' b))
>                                
>    td  _ ()  = (0,0) -- only core alignment counts
>    ti () _   = (0,0)
>    h []      = []
>    h l       = [minimum l]


length algebra, counts the flowgram core alignment length. 3' teminal gaps not included.

> length_alg :: FG_Algebra Float Int
> length_alg = (nil, d, i, r, td, ti, h) where
>    nil          = 0
>    d _ _ _ _ s  = s + 4
>    i s _ _ _ _  = s + 4
>    r a s b      = s + 1
>    td  (i,j) () = 0     -- only core residues count towards alignment length
>    ti () (i,j)  = 0
>    h            = id

seq_length algebra:
counts the length of the implicit sequence alignment

> seq_length_alg :: FG_Algebra Float Int
> seq_length_alg = (nil, d, i, r, td, ti, h) where
>    nil          = 0
>    d a b c d s  = s + round' a + round' b + round' c + round' d
>    i s a b c d  = s + round' a + round' b + round' c + round' d
>    r a s b      = s + max (round' a) (round' b)
>    td  (i,j) () = 0     -- only core residues count towards alignment length
>    ti () (i,j)  = 0
>    h            = id


Algebra cross product:
Note: uses take 1 in choice function. Needs to be remove if first alg is not
      single value optimization.

> infixl ***
> (***) :: Eq answer1 =>
>          FG_Algebra alphabet answer1 ->
>          FG_Algebra alphabet answer2 ->
>          FG_Algebra alphabet (answer1, answer2)
> alg1 *** alg2 = (nil, d, i, r, td, ti, h) where
>    (nil1, d1, i1, r1, td1, ti1, h1) = alg1
>    (nil2, d2, i2, r2, td2, ti2, h2) = alg2
> 
>    nil            = (nil1, nil2)
>    d  a b c d (s1,s2)   = (d1 a b c d s1, d2 a b c d s2)
>    i    (s1,s2) a b c d = (i1 s1 a b c d, i2 s2 a b c d)
>    r  x (s1,s2) y       = (r1 x s1 y, r2 x s2 y)
>    td r ()              = (td1 r (), td2 r ())
>    ti () r              = (ti1 () r, ti2 () r)
>    h xs = take 1 [(x1,x2)| x1 <- nub $ h1 [ y1 | (y1,y2) <- xs],
>                            x2 <-       h2 [ y2 | (y1,y2) <- xs, y1 == x1]]

The yield grammar:

> flowgram_ali:: Bool -> FG_Algebra Float b -> Flowgram -> Flowgram -> [b]
> flowgram_ali banded alg inpX inpY = axiom alignment where
>   (nil, d, i, r, td, ti, h) = alg
> 
>   alignment = tabulated (
>               nil ><< empty                                                                                     |||
>               r   <<<                               xbase -~~ alignment ~~- ybase                               |||
>               d   <<< xbase ~~- xbase ~~- xbase ~~- xbase ~~! alignment                                         |||
>               i   <<<                                         alignment ~~- ybase ~~- ybase ~~- ybase ~~- ybase |||
>               ti  <<<            empty .~~ region                                                               |||
>               td  <<< region ~~. empty                                                                          ... h)

>               where infixl 7 ~~!
>                     (~~!) = (~~*) (4,4) 0

Bind input:

>   z          = mk (toInput inpX inpY)
>   xlen       = length inpX
>   ylen       = length inpY

>   (_,n)      = bounds z
>   tabulated' = table' xlen ylen  -- standard rectangular tabulation
>   tabulated
>       | banded == True = banded_table' band_width  xlen ylen -- use this table for banded DP 
>       | otherwise      = table' xlen ylen  -- standard rectangular tabulation
>	where band_width = min 20 (max 10 (xlen `div` 20)) -- 10 < width < 20	
>   xbase      = achar' z
>   ybase      = achar' z
>   region     = astring
>   axiom      = axiom' n

>   empty (i,j) =  [() | i == xlen && j == xlen] -- empty marks the end of the recursion

> check_args:: [String] -> IO ()
> check_args args = case length args of 
>              3          -> return ()
>              otherwise  -> putStr help_string
>                            >> exitWith ExitSuccess

> dispatch:: String -> FlowSignalDistrib -> [Flowgram] -> String
> dispatch mode error_profile flows = case mode of
>        "-align"           -> concat $ map (\a -> show a ++"\n") (map (align error_profile (head flows)) (tail flows))
>        "-score"           -> concat $ map (\a -> show a ++"\n") (map (score_it error_profile (head flows)) (tail flows)) 
>        "-flow-lengths"    -> show   $ map length (map (flow_to_seq def_fo) flows)
>        "-relscore"        -> concat $ map (\a -> show a ++"\n") (filter_with_flow error_profile flows)
>        "-relscore_pairid" -> concat $ map (\(a,b)-> show a ++"\t"++show b++"\n") (filter_with_flow' error_profile flows)
>        "-self"            -> concat $ map (\a -> show a ++"\n") (map (self_score_fast error_profile) flows)
>        "-gapless"         -> concat $ map (\a -> show a ++"\n") (map (no_align_score error_profile (head flows)) (tail flows))
>        "-h"               -> help_string
>        "-help"            -> help_string
>        otherwise          -> help_string

> main = do
>        args <- getArgs
>        check_args args
>        let [mode, datfile, infile1] = args                                        
>        s             <- readFile infile1
>        error_profile <- readFlowSignalDistrib 1000 datfile
>        let (desc:seqs) = lines s    
>     --       flows = parseSFFfile (lines s)     -- .sff.txt format
>            flows = map (map read) (map (drop 2) (map words seqs)) ::[Flowgram]  --pyronoise format
>            result = dispatch mode error_profile flows   
>	 putStr (result++ "\n")

> help_string = unlines ["Usage: FlowgramAli_4frame mode error_profile input_file",
>                "where mode:",
>                "\t-align   \t\t Align all flowgrams in input against first flowgram in input",
>                "\t-score   \t\t Only compute alignment score",
>                "\t-flow-lengths\t\t print translated nucleotide length if flowgrams in input" ,
>                "\t-relscore\t\t Compute lenghth normalized alignment score",
>                "\t-relscore_pairid\t Compute length normalized alignment score and report %pair id of aligned seqs",
>                "\t-self    \t\t Fast self-alignment score for all flowgrams in input",
>                "\t-gapless\t\t Fast gapless alignment of first flowgram in input against rest in input"]

-----------------------------------
Helper functions

> format::(Int,(String,String)) -> String
> format (a,(b,c)) = "Score: " ++ show a ++ "\n" ++ b ++ "\n" ++ c

> toInput:: [a] -> [a] -> [a]
> toInput x y = x ++ (reverse y)

> align arr x y       = flowgram_ali True ((score arr) *** (prettyprint (mk (toInput x y)))) x y
> align' arr x y      = flowgram_ali True ((score arr) *** (prettyprint' arr (mk (toInput x y)))) x y

> score_it arr        = flowgram_ali  True (score arr)
> score_rel arr x y   = s / fromIntegral l
>     where [(s,l)] = flowgram_ali  True (score arr *** length_alg) x y
> score_rel_mismatch arr x y   = (s / fromIntegral l, 1 - fromIntegral m / fromIntegral seq_len)
>     where [(s,(l,(m,seq_len)))] = flowgram_ali  True (score arr *** (length_alg *** mismatch_seqlen)) x y

 compute alignment score if flowgram against itself using explicit alignment

> self_score_slow:: FlowSignalDistrib -> Flowgram -> Float
> self_score_slow arr x = case score_it arr x x of 
>                     [x] -> x
>                     otherwise -> error ("Bad self_score")

compute alignment score by doing fast, gapless self alignment

> self_score_fast :: FlowSignalDistrib -> Flowgram -> Float
> self_score_fast arr [] = 0
> self_score_fast arr (x:xs) = arr!(round' x, floor( x * 100)) + self_score_fast arr xs

compute alignments score of two flowgrams using gapless alignment

> no_align_score:: FlowSignalDistrib -> Flowgram -> Flowgram -> Float
> no_align_score arr a b = nas a b / fromIntegral(min (length a) (length b))
>      where
>      nas::  Flowgram -> Flowgram -> Float
>      nas [] [] = 0
>      nas (x:xs) (y:ys) = arr!(round' x, floor( y * 100)) + nas xs ys
>      nas  _ _ = 0

> filter_with_flow:: FlowSignalDistrib -> [Flowgram] -> [Float]
> filter_with_flow arr [] = []
> filter_with_flow arr (a:as) =  map (score_rel arr a) as

> filter_with_flow':: FlowSignalDistrib -> [Flowgram] -> [(Float,Float)]
> filter_with_flow' arr [] = []
> filter_with_flow' arr (a:as) =  map (score_rel_mismatch arr a) as
