%==========================================================================
%
% batch_hatch  Runs a hatch filter for a batch of GNSS data.
%
%   srho = hatch_filter(rho,phi)
%
% Author: Tamas Kis
% Last Update: 2022-03-07
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   rho     - (n×T double) pseudorange measurements [m]
%   phi     - (n×T double) carrier phase measurements [m]
%   SVID    - (n×T double) SVIDs of GPS satellites corresponding to 
%             measurements [m]
%
% -------
% OUTPUT:
% -------
%   rho     - (n×T double) smoothed carrier phase measurements [m]
%
% -----
% NOTE:
% -----
%   --> n = number of pseudorange measurements at each sample time
%   --> T = number of sample times
%
%==========================================================================
function rho = batch_hatch(rho,phi,SVID)

    % number of pseudorange measurements at each sample time
    n = size(rho,1);

    % number of sample times
    T = size(rho,2);

    % keeps track of hatch iterations
    M = zeros(n,1);

    % loop over all sample times
    for k = 2:T

        % loop over all GNSS satellites
        for i = 1:n
            
            % SVID of current satellite
            SVID_cur = SVID(i,k);
            
            % finds data index corresponding to satellite with same
            % SVID at previous iteration
            i_prev = find(SVID_cur == SVID(:,k-1),1);
            
            % determines if current satellite is a "new" satellite
            new_sat = isempty(i_prev);
            
            % initializes hatch filter if new satellite by resetting
            % counter
            if new_sat || (k == 2)
                M(i) = 0;

            % smooths pseudorange
            else
                rho(i,k) = hatch_filter(rho(i,k),rho(i_prev,k-1),...
                    phi(i,k),phi(i_prev,k-1),M(i));
            end

            % increments counter for hatch filter
            M = M+1;
            
        end

    end

end