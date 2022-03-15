%==========================================================================
%
% k2t  Time from sample number.
%
%   tk = k2t(k,dt)
%   tk = k2t(k,dt,t0)
%
% See also t2k.
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
%   k       - (1×1 double) current sample number
%   dt      - (1×1 double) time step
%   t0      - (1×1 double) (OPTIONAL) initial time (defaults to 0)
%
% -------
% OUTPUT:
% -------
%   tk      - (1×1 double) current sample time
%
%==========================================================================
function tk = k2t(k,dt,t0)
    
    % defaults initial time to 0 if not specified
    if (nargin < 3) || isempty(t0)
        t0 = 0;
    end
    
    % current sample time
    tk = t0+k*dt;
    
end