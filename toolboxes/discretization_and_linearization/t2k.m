%==========================================================================
%
% t2k  Sample number from time.
%
%   k = k2t(t,dt)
%   k = k2t(t,dt,t0)
%
% See also k2t.
%
% Author: Tamas Kis
% Last Update: 2022-02-25
%
% REFERENCES:
%   [1] TODO
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   t       - (1×1 double) current time
%   dt      - (1×1 double) time step
%   t0      - (1×1 double) (OPTIONAL) initial time (defaults to 0)
%
% -------
% OUTPUT:
% -------
%   k       - (1×1 double) current sample number
%
%==========================================================================
function k = t2k(t,dt,t0)
    
    % defaults initial time to 0 if not specified
    if (nargin < 3) || isempty(t0)
        t0 = 0;
    end
    
    % current sample number
    k = (t-t0)/dt;
    
end