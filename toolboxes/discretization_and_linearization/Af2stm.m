%==========================================================================
%
% Af2stm  State transition matrix from continuous dynamics Jacobian and
% dynamics equation.
%
%   [Phi,xk1] = Af2stm(A,f,xk,uk,tk,dt)
%   [Phi,xk1] = Af2stm(A,f,xk,[],[],dt)
%
% See also f2stm_num.
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
%   A       - (1×1 function_handle) continuous dynamics Jacobian, A(x,u,t)
%             (A : ℝⁿ×ℝᵐ×ℝ → ℝⁿˣⁿ)
%   f       - (1×1 function_handle) continuous dynamics equation,
%             dx/dt = f(x,u,t) (f : ℝⁿ×ℝᵐ×ℝ → ℝⁿ)
%   xk      - (n×1 double) state vector at kth sample time, xₖ
%   uk      - (m×1 double) (OPTIONAL) control input at kth sample time, uₖ
%   tk      - (1×1 double) (OPTIONAL) time at kth sample time, tₖ
%   dt      - (1×1 double) time step
%
% -------
% OUTPUT:
% -------
%   Phi     - (n×n double) state transition matrix, Φ(tₖ₊₁,tₖ)
%   xk1     - (n×1 double) state vector at (k+1)th sample time, Φ(tₖ₊₁,tₖ)
%
% -----
% NOTE:
% -----
%   --> If you would not like to specify "uk" or "tk", you can input them
%       as the empty vector "[]".
%
%==========================================================================
function [Phi,xk1] = Af2stm(A,f,xk,uk,tk,dt)
    
    % state dimension
    n = length(xk);
    
    % initial condition for augmented STM (at time tₖ)
    Psik = [eye(n),xk];
    
    % function handle for dΨ/dt = [AΦ f]
    dPsidt = @(t,Psi) [A(Psi(:,n+1),uk,t)*Psi(:,1:n),f(Psi(:,n+1),uk,t)];
    
    % solve for Ψ(tₖ₊₁,tₖ)
    Psi = RK4_step(dPsidt,tk,Psik,dt);
    
    % extract Φ(tₖ₊₁,tₖ) and xₖ₊₁ from the augmented STM Ψ(tₖ₊₁,tₖ)
    Phi = Psi(:,1:n);
    xk1 = Psi(:,n+1);
    
end