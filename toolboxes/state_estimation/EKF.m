%==========================================================================
%
% EKF  Extended Kalman filter.
%
%   [xk,Pk,z_pre,z_post,F_prev,Hk] = EKF(x_prev,P_prev,u_prev,yk,k,fd,...
%       hd,F,H,Q_prev,Rk)
%
% Author: Tamas Kis
% Last Update: 2022-03-05
%
% REFERENCES:
%   [1] TODO
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x_prev  - (n×1 double) state estimate at previous sample time
%   P_prev  - (n×n double) error covariance at previous sample time
%   u_prev  - (m×1 double) control input at previous sample time
%   yk      - (p×1 double) measurement at current sample time
%   k       - (1×1 double) current sample number
%   fd      - (1×1 function_handle) discrete nonlinear dynamics equation,
%             xₖ₊₁ = fd(xₖ,uₖ,k) (fd : ℝⁿ×ℝᵐ×ℤ → ℝⁿ)
%   hd      - (1×1 function_handle) discrete nonlinear measurement 
%             equation, yₖ = hd(xₖ,k) (fd : ℝⁿ×ℤ → ℝᵖ)
%   F       - (1×1 function_handle) Fₖ = F(xₖ,uₖ,k) --> discrete dynamics
%             Jacobian (F : ℝⁿ×ℝᵐ×ℤ → ℝⁿˣⁿ)
%   H       - (1×1 function_handle) Hₖ = H(xₖ,uₖ) --> discrete measurement
%             Jacobian (H : ℝⁿ×ℤ → ℝᵖˣⁿ)
%   Q_prev  - (n×n double) process noise covariance at previous sample time
%   Rk      - (p×p double) measurement noise covariance at current sample
%             time
%
% -------
% OUTPUT:
% -------
%   xk      - (n×1 double) a posteriori state estimate
%   Pk      - (n×n double) a posteriori error covariance
%   z_pre   - (p×1 double) pre-fit measurement residual
%   z_post  - (p×1 double) post-fit measurement residual
%   F_prev  - (n×n double) discrete dynamics Jacobian at previous sample 
%             time
%   Hk      - (p×m double) discrete measurement Jacobian at current sample
%             time
%
%==========================================================================
function [xk,Pk,z_pre,z_post,F_prev,Hk] = EKF(x_prev,P_prev,u_prev,yk,k,...
    fd,hd,F,H,Q_prev,Rk)
    
    % predict step (time update)
    [x_pred,P_pred,F_prev] = EKF_predict(x_prev,P_prev,u_prev,k,fd,F,...
        Q_prev);

    % update step (measurement update)
    [xk,Pk,z_pre,z_post,Hk] = EKF_update(x_pred,P_pred,yk,k,hd,H,Rk);
    
end