import Control.Applicative
import Data.List
import Graphics.Rendering.Chart.Simple
import Data.Random.Normal

import TupleKalman
import TupleAlgebra

-- Model.
f = ((1,1),
     (0,1)) :: F    -- Straight line.
h = (1,0) :: H       -- Direct measurement.
-- Noise guess.
q = ((1e-5,0),
     (0,1e-5))  -- Process noise covariance
r = 0.5^2       -- Measurement noise covariance.
-- State guess.
x0 = (0,0)
p0 = ((0.1^2,0),
      (0,0.1^2))

-- Observations
x0_true = (5,-0.1)

vs = mkNormals 1940141
zs = snd $ mapAccumL (observation f h) x0_true vs

observation :: F -> H -> X -> R -> (X,Z)
observation f h x_ v = (x', observe_z h x' + v)
  where
    x' = fst $ predict f undefined (x_,undefined)
    observe_z h x = h `dot` x

(xs,ps) = unzip $ discreteKF f q h r (x0,p0) zs

main = do
  plotWindow [1..100::Double]
             zs "o"
             (tail $ fst <$> fst <$> discreteKF f q h r (x0,p0) zs)
             (tail $ zipWith (\x p -> fst x + sqrt (fst $ fst p)) xs ps)
             vs "o"  -- yellow shit
             (tail $ zipWith (\x p -> fst x - sqrt (fst $ fst p)) xs ps)
