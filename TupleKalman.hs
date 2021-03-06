module TupleKalman where

import TupleAlgebra

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


type X = (Double,Double)
type K = X
type H = X
type P = ((Double,Double),
          (Double,Double))
type F = P
type Q = P
type R = Z
type T = Z
type Z = Double
type Y = Z
type S = Z

-- Prediction
predict_x :: F -> X -> X
predict_x f x_ = f .*^ x_
predict_p :: F -> Q -> P -> P
predict_p f q p_ = (f .*. p_ .*. tr f) .+. q
-- | Predict state and covariance one step.
predict :: F -> Q -> (X,P) -> (X,P)
predict f q (x_,p_) = (predict_x f x_, predict_p f q p_)

-- Observation
innovation_y :: H -> X -> Z -> Y
innovation_y h x' z = z - (h`dot`x')
innovationCovariance_s :: H -> R -> P -> S
innovationCovariance_s h r p' = (h ^*. p' `dot` h) + r
gain_k :: H -> P -> S -> K
gain_k h p' s = (p' .*^ h) ^/ s  -- Optimal Kalman gain.

-- Update
update_x :: K -> Y -> X -> X
update_x k y x' = x' ^+^ (k ^* y)
update_p :: H -> K -> P -> P
update_p h k p' = (i .-. (k `dyad` h)) .*. p'
-- | Update predicted state and covariance using the given gain and innovation.
update' :: H -> K -> Y -> (X,P) -> (X,P)
update' h k y (x',p') = (update_x k y x', update_p h k p')
-- | Update predicted state and covariance using the given observation.
update :: H -> R -> Z -> (X,P) -> (X,P)
update h r z (x',p') = (update_x k y x', update_p h k p')
  where
    y = innovation_y h x' z
    s = innovationCovariance_s h r p'
    k = gain_k h p' s

-- | Predicts state and covariance and then updated the predicted state
-- and covariance using the given observation.
predictUpdate :: F -> Q -> H -> R -> (X,P) -> Z -> (X,P)
predictUpdate f q h r (x_,p_) z = update h r z $ predict f q (x_,p_)

-- | Discrete Kalman filter where the state transition and observation
-- matrices are constant (do not change between time steps).
discreteKF :: F -> Q -> H -> R -> (X,P) -> [Z] -> [(X,P)]
discreteKF f q h r = scanl (predictUpdate f q h r)

-- | Discrete Kalman filter where the state transition matrix is a
-- function of the time that has passed since the old state (i.e.
-- the duration of prediction). The observations must be paired with
-- the time difference since the preceeding observation.
deltaTimeKF :: (T->F) -> Q -> H -> R -> (X,P) -> [(T,Z)] -> [(X,P)]
deltaTimeKF f q h r  = scanl (\xp (dt,z) -> predictUpdate (f dt) q h r xp z)
