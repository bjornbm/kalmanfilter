{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables #-}

{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE UndecidableInstances #-}

module KalmanStatic where

import Numeric.Units.Dimensional.LinearAlgebra
import Numeric.Units.Dimensional.Prelude
import Data.HList
import qualified Prelude

-- Naming:
--
--   x  state
--   p  state covariance
--
--   f  state transition model
--   w  process noise
--   q  process noise covariance
--
--   b  control input model
--   u  control input vector
--
--   h  observation model
--   z  measurement
--   v  measurement noise
--   r  measurement noise covariance
--
--   y  innovation (or residual)
--   s  innovation covariance
--   k  optimal kalman gain
--
-- Convention:
--
--   x_  old (x-minus)
--   x'  predicted (x-pr)
--   x   updated


-- Convenience types.
type D d = Quantity d Double
type V v = Vec      v Double
type M m = Mat      m Double


-- Prediction
-- ----------
--  F -> X_ -> X'
predict_x' :: MatrixVector f x x => M f -> V x -> V x
predict_x' f x_ = f |*< x_

--  F -> Q -> P_ -> P'
predict_p' :: (Transpose f f', MatrixMatrix f p m, MatrixMatrix m f' p)
           => M f -> M p -> M p -> M p
predict_p' f q p_ = f |*| p_ |*| transpose f |+| q

-- | Predict state and covariance one step.
--  F -> Q -> (X_,P_) -> (X',P')
predict :: ( MatrixVector f x x, Transpose f f'
           ,  MatrixMatrix f p m, MatrixMatrix m f' p )
        => M f -> M p -> (V x, M p) -> (V x, M p)
predict f q (x_,p_) = (predict_x' f x_, predict_p' f q p_)


-- Observation
-- -----------
--  H -> X' -> Z -> Y
innovation_y :: DotProduct h x z => V h -> V x -> D z -> D z
innovation_y h x' z = z - (h >.< x')

--  H -> R -> P' -> S
innovationCovariance_s h r p' = h >*| p' >.< h + r

--  H -> S -> P' -> K
gain_k h s p' = p' |*< h >/ s  -- Optimal Kalman gain.


-- Update
-- ------
--  K -> X' -> Y -> X
update_x :: (HMap (MulD,y) k x) => V k -> V x -> D y -> V x
update_x k x' y = x' >+< (k >* y)

--  H -> K -> P' -> P
update_p h k p' = (identity |-| k `dyadicProduct` h) |*| p'

-- | Update predicted state and covariance using the given gain and
-- innovation. This is presented here for completeness. Normally you
-- will want to use 'update' instead.
--  H -> K -> (X',P') -> y -> (X,P)
update0 h k (x',p') y = (update_x k x' y, update_p h k p')

-- | Update predicted state and covariance using the given observation.
--  H -> R -> (X',P') -> Z -> (X,P)
update h r (x',p') z = (update_x k x' y, update_p h k p')
  where
    y = innovation_y h x' z
    s = innovationCovariance_s h r p'
    k = gain_k h s p'


-- Combined prediction and update
-- ------------------------------
-- | Predicts state and covariance and then updated the predicted state
-- and covariance using the given observation.
--  F -> Q -> H -> R -> (X_,P_) -> Z -> (X,P)
predictUpdate f q h r (x_,p_) z = update h r (predict f q (x_,p_)) z


-- Kalman filters
-- --------------
-- | Discrete Kalman filter where the state transition and observation
-- matrices are constant (do not change between time steps).
--  F -> Q -> H -> R -> (X0,P0) -> [Z] -> [(X,P)]
discreteKF f q h r = scanl (predictUpdate f q h r)

-- | Discrete Kalman filter where the state transition matrix is a
-- function of the time that has passed since the old state (i.e.
-- the duration of prediction). The observations must be paired with
-- the time difference since the preceeding observation.
--  (T->F) -> (T->Q) -> H -> R -> (X0,P0) -> [(T,Z)] -> [(X,P)]
deltaTimeKF f q h r = scanl (\xp (dt,z) -> predictUpdate (f dt) (q dt) h r xp z)
