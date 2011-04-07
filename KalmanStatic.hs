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
type V v  = Vec v  Double
type M vs = Mat vs Double


-- Prediction
--predict_x :: MatrixVector vs v v => F vs -> X v -> X v
predict_x :: MatrixVector f x x => M f -> V x -> V x
predict_x f x_ = f |*< x_
--predict_p :: F -> Q -> P -> P
predict_p :: (Transpose f f', MatrixMatrix f p m, MatrixMatrix m f' p) => M f -> M p -> M p -> M p
predict_p f q p_ = f |*| p_ |*| transpose f |+| q
-- | Predict state and covariance one step.
--predict :: F -> Q -> (X,P) -> (X,P)
predict :: (MatrixVector f x x, Transpose f f', MatrixMatrix f p m, MatrixMatrix m f' p)
        => M f -> M p -> (V x, M p) -> (V x, M p)
predict f q (x_,p_) = (predict_x f x_, predict_p f q p_)

-- Observation
--innovation_y :: H -> X -> Z -> Y
innovation_y h x' z = z - (h >.< x')
--innovationCovariance_s :: H -> R -> P -> S
innovationCovariance_s h r p' = h >*| p' >.< h + r
--gain_k :: H -> P -> S -> K
gain_k h p' s = p' |*< h >/ s  -- Optimal Kalman gain.

-- Update
--update_x :: K -> Y -> X -> X
update_x k y x' = x' >+< (k >* y)
--update_p :: H -> K -> P -> P
update_p h k p' = (identity |-| k `dyadicProduct` h) |*| p'

-- | Update predicted state and covariance using the given gain and innovation.
--update' :: H -> K -> Y -> (X,P) -> (X,P)
update' h k y (x',p') = (update_x k y x', update_p h k p')
-- | Update predicted state and covariance using the given observation.
--update :: H -> R -> Z -> (X,P) -> (X,P)
update h r z (x',p') = (update_x k y x', update_p h k p')
  where
    y = innovation_y h x' z
    s = innovationCovariance_s h r p'
    k = gain_k h p' s

-- | Predicts state and covariance and then updated the predicted state
-- and covariance using the given observation.
--predictUpdate :: F -> Q -> H -> R -> (X,P) -> Z -> (X,P)
predictUpdate f q h r (x_,p_) z = update h r z $ predict f q (x_,p_)

-- | Discrete Kalman filter where the state transition and observation
-- matrices are constant (do not change between time steps).
--discreteKF :: F -> Q -> H -> R -> (X,P) -> [Z] -> [(X,P)]
discreteKF f q h r = scanl (predictUpdate f q h r)

-- | Discrete Kalman filter where the state transition matrix is a
-- function of the time that has passed since the old state (i.e.
-- the duration of prediction). The observations must be paired with
-- the time difference since the preceeding observation.
--deltaTimeKF :: (T->F) -> Q -> H -> R -> (X,P) -> [(T,Z)] -> [(X,P)]
deltaTimeKF f q h r  = scanl (\xp (dt,z) -> predictUpdate (f dt) q h r xp z)
