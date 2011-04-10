{-# LANGUAGE TypeOperators #-}


import Control.Applicative
import Data.List hiding (sum)
import Graphics.Rendering.Chart.Simple
import Data.Random.Normal (mkNormals)

import KalmanStatic
import Numeric.Units.Dimensional.Prelude
import Numeric.Units.Dimensional.LinearAlgebra
import Numeric.NumType (Zero, Pos1, Pos2, Neg1, Neg2, Neg3, Neg4)
import qualified Prelude

-- Dimensionality of matrix and vector elements.
type DVelocityInv    = Dim Neg1 Zero Pos1 Zero Zero Zero Zero
type DVelocitySq     = Dim Pos2 Zero Neg2 Zero Zero Zero Zero
type DAccelerationSq = Dim Pos2 Zero Neg4 Zero Zero Zero Zero
type DLengthSq       = Dim Pos2 Zero Zero Zero Zero Zero Zero
type DLengthVel      = Dim Pos2 Zero Neg1 Zero Zero Zero Zero
type DTimeSq         = Dim Zero Zero Pos2 Zero Zero Zero Zero
type DTimeInv        = Dim Zero Zero Neg1 Zero Zero Zero Zero

-- Be specific with types to help the type checker.
type X = Vec  (DLength   :*.DVelocity)   Double
type P = Mat ((DLengthSq :*.DLengthVel )
          :*. (DLengthVel:*.DVelocitySq)) Double
type F = Mat ((DOne    :*.DTime)
          :*. (DTimeInv:*.DOne )) Double
type Q = P

type H = Vec (DOne:*.DTime) Double
type Z = Length Double
type R = Quantity DLengthSq Double

type Y = Z
type S = Quantity DLengthSq Double
type K = X

type T = Time Double

type G = Vec (DTimeSq:*.DTime) Double
type A = Acceleration Double


-- Kalman filter setup (model)
-- ===========================

-- | State transition model.
f :: T -> F
f dt = (_1             <:. dt) |:.
       (0*~second^neg1 <:. _1)

-- | Observation model.
h :: H
h = _1 <:. 0*~second  -- Direct measurement.

-- | Process noise covariance guess (we assume sigma_a is accurately known).
q :: T -> Q
q dt = ((dt^pos4     <:. dt^pos3 /_2)  |:.
        (dt^pos3 /_2 <:. dt^pos2    )) |* sigma_a^pos2

-- | Measurement noise covariance guess (we assume sigma_z is accurately known).
r :: R
r = sigma_z^pos2

-- | Initial state guess (we assume perfect knowledge).
x0 :: X
x0 = x0_true
p0 :: P
p0 = (0*~meter^pos2          <:. 0*~(meter^pos2/second)     ) |:.
     (0*~(meter^pos2/second) <:. 0*~(meter^pos2/second^pos2))


-- Physical reality
-- ================

-- | Initial state.
x0_true :: X
x0_true = 0*~meter <:. 0*~(meter/second)

-- | Unmodeled acceleration impact
g :: T -> G
g dt = dt^pos2 /_2 <:. dt

-- | Propagation of the true state for the given time step an mean
-- acceleration.
propagate_x :: X -> (T,A) -> X
propagate_x x_true_ (dt,a) = predict_x (f dt) x_true_ >+< g dt >* a

-- | Variance of mean acceleration.
sigma_a :: Quantity DAcceleration Double
sigma_a = 0.1*~(meter/second^pos2)

-- | Mean accelerations.
as :: [A]
as = (* sigma_a) <$> mkNormals 1847101 *~~ one


-- Sampling
-- ========

-- | Assume we sample once per second.
dts :: [T]
dts = repeat 1 *~~ second

-- | Sampling noise standard deviation.
sigma_z :: Z
sigma_z = 1.0 *~ meter

-- | The measurement errors.
vs :: [Z]
vs = (* sigma_z) <$> mkNormals 908714 *~~ one

-- | True state at t0, t1...
xs_true :: [X]
xs_true = scanl propagate_x x0_true (zip dts as)

-- | Measurements at t1, t2... (note how we use only the tail of xs_true).
zs :: [Z]
zs = zipWith (+) (map (h >.<) $ tail xs_true) vs


-- Filtering
-- =========
xps :: [(X,P)]
xps = deltaTimeKF f q h r (x0,p0) (zip dts zs)

xs :: [X]
xs = fst <$> xps

ps :: [P]
ps = snd <$> xps


main = do
  plotWindow ([1..120] :: [Double])
             (map (>!!zero) xs_true /~~ meter)
             (0:zs/~~meter) "o"
             (map (>!!pos1) xs_true /~~ (deci meter/second))
             (repeat 1 :: [Double])  -- yellow
             (as /~~ (deci meter/second^pos2))
             (map (>!!zero) xs /~~ meter)
             (map (>!!pos1) xs /~~ (deci meter/second))
             (map (sqrt . (>!!zero) . rowHead) ps /~~ meter)
             (map (sqrt . (>!!pos1) . rowHead . rowTail) ps /~~ (meter/second))
