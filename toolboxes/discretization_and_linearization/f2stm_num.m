%==========================================================================
%
% f2stm_num  State transition matrix from continuous dynamics equation.
%
%   Phi = f2stm_num(f,xk,uk,tk,dt)
%   Phi = f2stm_num(f,xk,[],[],dt)
%
% See also A2stm_num.
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
%
% -----
% NOTE:
% -----
%   --> If you would not like to specify "uk" or "tk", you can input them
%       as the empty vector "[]".
%
%==========================================================================
function Phi = f2stm_num(f,xk,uk,tk,dt)
    
    % function handle for dynamics Jacobian
    A = @(x,u,t) f2A_num(f,x,u,t);

    % solve for Φ(tₖ₊₁,tₖ)
    Phi = Af2stm(A,f,xk,uk,tk,dt);
    
end