%==========================================================================
%
% imod "Complexified" version of the modulo operator.
%
%   r = imod(a,n)
%
% See also mod.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-03-09
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TOOLBOX DOCUMENTATION:
% https://tamaskis.github.io/Numerical_Differentiation_Toolbox-MATLAB/
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/documentation/Numerical_Differentiation_using_the_Complex_Step_Approximation.pdf
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x       - (n×1 complex) input argument
%   y       - (n×1 complex) input argument
%
% -------
% OUTPUT:
% -------
%   z       - (1×1 complex) dot product of x and y
%
%==========================================================================
function r = imod(a,n)
    r = a-floor(a/n)*n;
end