{-# LANGUAGE TypeOperators #-}


import Control.Applicative
import Data.List
import Graphics.Rendering.Chart.Simple
import System.Random

import Vector
import Matrix
import KalmanStatic
import Numeric.Units.Dimensional.Prelude
import qualified Prelude

-- Be specific with types to help the type checker.
type X = Vec (DOne :*. DOne)  Double
type P = Mat ((DOne :*. DOne) :*.
              (DOne :*. DOne)) Double
type F = P
type Q = P

type H = X
type Z = Quantity DOne Double
type R = Z

type Y = Z
type S = Z
type K = X

type T = Z

-- Helper for easier matrix construction.
fromTuples (v1,v2) = consRow   (fromTuple v1) $
                     rowMatrix (fromTuple v2)

-- Model.
f = fromTuples ((_1,_1),
                (_0,_1)) :: F    -- Straight line.
h = fromTuple (_1,_0) :: H       -- Direct measurement.
-- Noise guess.
q = fromTuples ((1e-5*~one,_0),
                (_0,1e-5*~one)) :: Q  -- Process noise covariance
r = (0.5*~one)^pos2       -- Measurement noise covariance.
-- State guess.
x0 = fromTuple (_0,_0) :: X
p0 = fromTuples (((0.1*~one)^pos2,_0),
                 (_0,(0.1*~one)^pos2)) :: P

-- Observations
x0_true = fromTuple (_5,(-0.1)*~one) :: X

ts = [0..100] *~~ one
vs = map cos ts
zs = snd $ mapAccumL (observation f h) x0_true vs

observation :: F -> H -> X -> R -> (X,Z)
observation f h x_ v = (x', observe_z h x' + v)
  where
    x' = fst $ predict f (undefined::Q) (x_,undefined::P)
    observe_z h x = h `dotProduct` x

(xs,ps) = unzip $ discreteKF f q h r (x0,p0) zs

main = do
  plotWindow (ts/~~one)
             (zs/~~one) "o"
             ((/~~one) $ vHead <$> fst <$> discreteKF f q h r (x0,p0) zs)
