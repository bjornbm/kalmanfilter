import Numeric.AD
import Graphics.Rendering.Chart.Simple
import System.Random
import Control.Applicative

--   z  measurement
--   x  state
--   p  state covariance
--   q  process variance (noise)
--   r  measurement variance (noise)

-- Convention:
--   x_  old (x-minus)
--   x'  predicted (x-pr)
--   x   updated

measurement_z x v = x + v

predict_x x_ = x_
predict_p q p_ = p_ + q

update_k r p'   = p' / (p' + r)
update_x k x' z = x' + k * (z - x')
update_p k p'   = (1 - k) * p'



type X = Double

-- Noise
vs :: [Double]
vs = map cos [0,1..100]

x_true = -2
measurements = measurement_z x_true <$> vs

r = 0.5^2
q = 1e-5

kfilter :: X -> X -> (X,X) -> X -> (X,X)
kfilter q r (x_,p_) z = (update_x k x' z, update_p k p')
  where
    x' = predict_x x_
    p' = predict_p q p_
    k  = update_k r p'

states q r xp0 zs = scanl (kfilter q r) xp0 zs

(xs,ps) = unzip $ states q r (0,5) measurements

main = plotWindow [0..100::Double]
       measurements "o"
       xs
       (zipWith (\x p -> x + sqrt p) xs ps)
       vs "o"  -- yellow shit
       (zipWith (\x p -> x - sqrt p) xs ps)
