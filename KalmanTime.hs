--import Numeric.AD
import Graphics.Rendering.Chart.Simple
import System.Random
import Control.Applicative

-- Naming:
--   z  measurement
--   x  state
--   p  state covariance
--   q  process variance (noise)
--   r  measurement variance (noise)
--
-- Convention:
--   x_  old (x-minus)
--   x'  predicted (x-pr)
--   x   updated

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
x -. ((a,b),(c,d)) = ((x-a,x-b),(x-c,x-d))
i = ((1,0),
     (0,1))

-- The model h(t)
model t = (1,t)  -- h(t)

measurement_z :: (T -> X) -> X -> T -> Z
measurement_z model x_true t = (model t `dot` x_true)

predict_x :: T -> X -> X
predict_x t x_ = x_
predict_p :: Q -> T -> P -> P
predict_p q t p_ = p_ .+. q

update_k :: H -> R -> P -> K
update_k h r p' = (p' .*^ h) ^/ ((h ^*. p' `dot` h) + r)
update_p :: H -> K -> P -> P
update_p h k p'   = (i .-. (k `dyad` h)) .*. p'
update_x :: (T -> H) -> T -> K -> X -> Z -> X
update_x model t k x' z =  k ^* (z - measurement_z model x' t)



type X = (Double,Double)
type K = X
type H = X
type P = ((Double,Double),
          (Double,Double))
type Q = P
type R = Z
type T = Z
type Z = Double

-- Times
ts :: [Double]
ts = [0..100]
-- Measurement noise
vs :: [Double]
vs = map cos [0..100]
x_true = (2,0.1)
-- Measurements
measurements = zipWith (+) (measurement_z model x_true <$> ts) vs
-- Initial guess
x0 = (0,0) :: X
p0 = ((1,1),
      (1,1)) :: P
r  = 0.5^2 :: R
q  = ((1e-5,    0),
      (   0, 1e-5)) :: Q


kfilter :: Q -> R -> (X,P) -> (T,Z) -> (X,P)
kfilter q r (x_,p_) (t,z) = (x,p)
  where
    x' = predict_x t x_
    p' = predict_p q t p_
    k  = update_k (model t) r p'
    x  = update_x model t k x' z
    p  = update_p (model t) k p'

states q r xp0 zs = scanl (kfilter q r) xp0 zs

(xs,ps) = unzip $ states q r (x0,p0) (zip ts measurements)

main = do
  print ps
  plotWindow ts
             measurements "o"
             (zipWith ($) (measurement_z model <$> xs) ts)
             --(zipWith (\x p -> x + sqrt p) xs ps)
             --vs "o"  -- yellow shit
             --(zipWith (\x p -> x - sqrt p) xs ps)
-- -}
