----------------------------------------------------------------------
ADP combinators and functions from:

R. Giegerich, C. Meyer and P. Steffen. Towards a discipline of dynamic
programming.

----------------------------------------------------------------------

> module ADPCombinators where
> import Data.Array

Lexical parsers
----------------

> type Subword  = (Int,Int)
> type Parser b = Subword -> [b]

> empty        :: Parser ()
> empty  (i,j) =  [() | i == j]

> acharSep'           ::  Array Int Char -> Char -> Parser Char
> acharSep' z s (i,j) =  [z!j | i+1 == j, z!j /= s] 

> achar'         :: Array Int a -> Parser a
> achar' z (i,j) = [z!j | i+1 == j]

> char'           ::  Eq a => Array Int a -> a -> Parser a
> char' z c (i,j) =  [c | i+1 == j, z!j == c]

> astring       :: Parser Subword
> astring (i,j) =  [(i,j) | i <= j]

> string'           :: Eq a => Array Int a -> [a] -> Parser Subword
> string' z s (i,j) = [(i,j)| and [z!(i+k) == s!!(k-1) | k <-[1..(j-i)]]]

Parser combinators
-------------------

> infixr 6 ||| 
> (|||)           :: Parser b -> Parser b -> Parser b
> (|||) r q (i,j) = r (i,j) ++ q (i,j)

> infix  8 <<<
> (<<<)           :: (b -> c) -> Parser b -> Parser c
> (<<<) f q (i,j) =  map f (q (i,j))

> infixl 7 ~~~
> (~~~)           :: Parser (b -> c) -> Parser b -> Parser c
> (~~~) r q (i,j) =  [f y | k <- [i..j], f <- r (i,k), y <- q (k,j)]  

> infix  5 ...
> (...)           :: Parser b -> ([b] -> [b]) -> Parser b
> (...) r h (i,j) = h (r (i,j))


> type Filter    =  (Int, Int) -> Bool
> with           :: Parser b -> Filter -> Parser b
> with q c (i,j) =  if c (i,j) then q (i,j)  else []

> axiom'        :: Int -> Parser b -> [b]
> axiom' l ax   =  ax (0,l) 

Tabulation
-----------

> -- two-dimensional tabulation
> table     :: Int -> Parser b -> Parser b
> table n q =  ((!) $ array ((0,0),(n,n))
>                    [((i,j),q (i,j)) | i<- [0..n], j<- [i..n]])

> -- two-dimensional tabulation with rectangular table, used e.g. in alignments
> -- instead of the large triangular matrix (n+m x n+m) as used in table, we
> -- store only the smaller (n x m) matrix and use a clever remapping.
> -- For n==m this reduces memory by four and even more for n/=m
> -- Note that his implementation requires no separating '$' symbol in between x and y
> table'     :: Int -> Int ->  Parser b -> Parser b
> table' n m q  = lookup $ array ((0,0),(n, m))
>                    [((i,j), q (i, n+j)) | i<- [0..n], j<- [0..m]]
>                    
>     where lookup:: Array (Int,Int) [b] -> Parser b
>           lookup arr (i,j) 
>               | i <= n && j >= n = arr!(i,j-n)
>               | otherwise = []

Banded tabulation for values around a band [(i,j) | j<-[i-k..i+k]]

> banded_table'     :: Int -> Int -> Int ->  Parser b -> Parser b
> banded_table' k n m q  = lookup $ array ((0,0),(n, m))
>                    [((i,j), q (i, n+j)) | i<- [0..n], j <- [max 0 (m-i-k) .. min (m-i+k) m]]
>                    
>     where lookup:: Array (Int,Int) [b] -> Parser b
>           lookup arr (i,j) 
>               | i <= n && j >= n && j > n+ (m-i)-k &&  j < n+(m-i) + k = arr!(i,j-n)
>               | otherwise = []


> -- one-dimensional tabulation; index j fixed
> listi :: Int -> Parser b -> Parser b
> listi n p = q $ array (0,n) [(i, p (i,n)) | i <- [0..n]] 
>    where
>    q t (i,j) = if j==n then t!i else []

> -- one-dimensional tabulation; index i fixed
> listj :: Int -> Parser b -> Parser b
> listj n p = q $ array (0,n) [(j, p (0,j)) | j <- [0..n]] 
>    where
>    q t (i,j) = if i==0 then t!j else []

> -- the most common listed type is listi (input read from left
> -- to right), so we define a default list here:
> list :: Int -> Parser b -> Parser b
> list = listi

Variants of the <<< and ~~~ Combinators
-------------------------------

> infix  8 ><<
> infixl 7 ~~, ~~*, *~~, *~*
> infixl 7 -~~, ~~-, +~~, ~~+, +~+, .~~, ~~.

The operator ><< is the special case of <<< for a nullary function f

> (><<)           :: c -> Parser b -> Parser c
> (><<) f q (i,j) =  [f|a <- (q (i,j))]

Subwords on left and right of an explicit length range.

> (~~) :: (Int,Int) -> (Int,Int) 
>      -> Parser (b -> c) -> Parser b -> Parser c
> (~~) (l,u) (l',u') r q (i,j) 
>      = [x y | k <- [max (i+l) (j-u') .. min (i+u) (j-l')],
>               x <- r (i,k), y <- q (k,j)]

Subwords of explicit length range and unbounded length on one or on
either side.

> (~~*) :: (Int,Int) -> Int 
>       -> Parser (a -> b) -> Parser a -> Parser b 
> (~~*) (l, u) l' r q (i, j) 
>       = [x y | k <- [(i + l) .. min (i + u) (j - l')], 
>                x <- r (i, k), y <- q (k, j)] 

> (*~~) :: Int -> (Int,Int) 
>       -> Parser (a -> b) -> Parser a -> Parser b 
> (*~~) l (l', u') r q (i, j) 
>       = [x y | k <- [max (i + l) (j - u') .. (j - l')], 
>           x <- r (i, k), y <- q (k, j)] 

> (*~*) :: Int -> Int 
>       -> Parser (a -> b) -> Parser a -> Parser b 
> (*~*) l l' r q (i, j) 
>       = [x y | k <- [(i + l) .. (j - l')], 
>                x <- r (i, k), y <- q (k, j)] 

Single character on the lefthand (respectively righthand) side

> (-~~)           :: Parser (b -> c) -> Parser b -> Parser c
> (-~~) q r (i,j) = [x y | i<j, x <- q (i,i+1), y <- r (i+1,j)]

> (~~-)           :: Parser (b -> c) -> Parser b -> Parser c
> (~~-) q r (i,j) = [x y | i<j, x <- q (i,j-1), y <- r (j-1,j)]

No character on the lefthand (respectively righthand) side

> (.~~)           :: Parser (b -> c) -> Parser b -> Parser c
> (.~~) q r (i,j) = [x y | i<j, x <- q (i,i), y <- r (i,j)]

> (~~.)           :: Parser (b -> c) -> Parser b -> Parser c
> (~~.) q r (i,j) = [x y | i<j, x <- q (i,j), y <- r (j,j)]

Nonempty sequence on the lefthand (respectively righthand) side

> (+~~)           :: Parser (b -> c) -> Parser b -> Parser c
> (+~~) r q (i,j) =  [f y | k <- [i+1..j], f <- r (i,k), y <- q (k,j)]  

> (~~+)           :: Parser (b -> c) -> Parser b -> Parser c
> (~~+) r q (i,j) =  [f y | k <- [i..j-1], f <- r (i,k), y <- q (k,j)]  

Nonempty sequence on either side

> (+~+)           :: Parser (b -> c) -> Parser b -> Parser c
> (+~+) r q (i,j) = [f y | k <- [(i+1)..(j-1)], f <- r (i,k), y <- q (k,j)]  


Create array from List
-----------------------

> mk :: [a] -> Array Int a
> mk xs = array (1,n) (zip [1..n] xs) where n = length xs


--------------------------------
Custom additions

> combiner low high = ((~~!), (~~!!))
>      where
>       (~~!) :: Parser (b -> c) -> Parser b -> Int -> (Int,Int) ->  [c]
>       (~~!) p q l (i,j)  = [x y |i+l < j, x <- p (i,i+l), y <- q (i+l,j)]
>       (~~!!) :: (Int -> (Int, Int) -> [(b -> c)]) -> Parser b -> Parser c
>       (~~!!) q r (i,j)  = [y z | l<-[low .. high], 0<= j-l, i+l < j,
>                            y <- (q l) (i,j-l), z <- r (j-l,j)]

