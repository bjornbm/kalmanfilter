module TupleAlgebra where

-- Algebra
(a,b)`dot`(c,d) = a*c + b*d
(a,b)`dyad`(c,d) = ((a*c,a*d),
                    (b*c,b*d))
(v1,
 v2) .*. ((a,b),
          (c,d)) = ((v1`dot`(a,c), v1`dot`(b,d)),
                    (v2`dot`(a,c), v2`dot`(b,d)))
v ^*. ((a,b),
       (c,d)) = (v`dot`(a,c), v`dot`(b,d))
(v1,v2) .*^ v = (v1`dot`v, v2`dot`v)
(a,b) ^* x = (a*x,b*x)
(a,b) ^/ c = (a/c,b/c)
((a,b),(c,d)) .+. ((e,f),(g,h)) = ((a+e,b+f),(c+g,d+h))
((a,b),(c,d)) .-. ((e,f),(g,h)) = ((a-e,b-f),(c-g,d-h))
(a,b) ^+^ (c,d) = (a+c,b+d)
x -. ((a,b),(c,d)) = ((x-a,x-b),(x-c,x-d))
i :: Num a => ((a,a),(a,a))
i = ((1,0),
     (0,1))
tr ((a,b),
    (c,d)) = ((a,c),
              (b,d))
