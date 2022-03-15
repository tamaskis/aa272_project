%==========================================================================
%
% measurement_jacobian_baseline  Continuous measurement Jacobian from 
% current state.
%
%   C = measurement_jacobian_baseline(x,t,prop,gps)
%   C = measurement_jacobian_baseline(x,t,prop,gps,c)
%
% Author: Tamas Kis
% Last Update: 2022-03-14
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x       - (16×1 double) state vector
%               -->   1-3. rc_eci - chief position resolved in ECI frame 
%                                   [m]
%               -->   4-6. vc_eci - chief inertial velocity resolved in ECI 
%                                   frame [m/s]
%               -->     7. bc     - chief clock bias [m]
%               -->     8. bc_dot - chief clock bias drift rate [m/s]
%               -->  9-11. rd_eci - deputy position resolved in ECI frame 
%                                   [m]
%               --> 11-14. vd_eci - deputy inertial velocity resolved in 
%                                   ECI frame [m/s]
%               -->    15. bd     - deputy clock bias [m]
%               -->    16. bd_dot - deputy clock bias drift rate [m/s]
%   t       - (1×1 double) simulation time [s]
%   prop    - (1×1 struct) propagator settings (see "initialize_propagator"
%             function for full definition)
%   gps     - (1×1 GPS_Constellation) GPS constellation object
%   l       - (1×1 double) (OPTIONAL) number of GNSS satellites used for 
%             measurements (defaults to 4)
%
% -------
% OUTPUT:
% -------
%   C       - (4c×16 double) continuous measurement Jacobian
%
%==========================================================================
function C = measurement_jacobian_baseline(x,t,prop,gps,l)
    
    % ---------------------
    % Unpacks state vector.
    % ---------------------
    
    % chief position [m] and inertial velocity [m/s] resolved in ECI frame
    rc_eci = x(1:3);
    vc_eci = x(4:6);
    
    % deputy position [m] and inertial velocity [m/s] resolved in ECI frame
    rd_eci = x(9:11);
    vd_eci = x(12:14);
    
    % -----------------------------
    % Timing and Earth orientation.
    % -----------------------------
    
    % time scales [MJD]
    [MJD_GPS,~,MJD_TT,MJD_UT1] = time_scales(t,prop.t0);
    
    % Earth orientation parameters for IAU2006/2000 CIO based theory
    [xp,yp,dX,dY,LOD] = eop_iau06(MJD_UT1,prop.data.eop);
    
    % rotation matrix (GCRF --> ITRF) and Earth angular velocity resolved
    % in the ITRF [rad/s] from IAU2006/2000 CIO based theory
    [~,R_eci2ecef,w_eci] = iau06(MJD_UT1,MJD_TT,xp,yp,dX,dY,...
        LOD,prop.data.XYs_iau06);
    
    % --------------------------------------------------------------------
    % Positions/velocities needed for computation of partials derivatives.
    % --------------------------------------------------------------------
    
    % chief position [m] resolved in ECEF frame
    rc_ecef = eci2ecef(rc_eci,vc_eci,w_eci,R_eci2ecef);
    
    % deputy position [m] resolved in ECEF frame
    rd_ecef = eci2ecef(rd_eci,vd_eci,w_eci,R_eci2ecef);

    % defaults "l" to 4 if not input
    if (nargin < 5) || isempty(l)
        l = 4;
    end
    
    % SVIDs of the "l" closest GPS satellites to chief
    SVID = gps.get_closest_SVIDs(rc_ecef,MJD_GPS,l);
    
    % ECEF positions [m] and velocities [m/s] of the "l" closest GPS 
    % satellites to chief
    r_sat_ecef = gps.get_ECEF_position(SVID,MJD_GPS);

    % ECI positions [m] and velocities [m/s] of the "l" closest GPS 
    % satellites to chief
    r_sat_eci = gps.get_ECI_position(SVID,MJD_GPS);
    v_sat_eci = gps.get_ECI_velocity(SVID,MJD_GPS);
    
    % -----------------------------------
    % Computation of partial derivatives.
    % -----------------------------------
    
    % partial derivatives of pseudorange for chief
    %[drhodr_c,drhodv_c] = dydX_pseudorange(rc_eci,r_sat_eci);
    [drhodr_c,drhodv_c] = dydX_pseudorange(rc_ecef,r_sat_ecef,...
        R_eci2ecef,w_eci);
    
    % partial derivatives of pseudorange rate for chief
    [drhodotdr_c,drhodotdv_c] = dydX_pseudorange_rate(rc_eci,r_sat_eci,...
        vc_eci,v_sat_eci);
    
    % partial derivatives of pseudorange for deputy
    [drhodr_d,drhodv_d] = dydX_pseudorange(rd_ecef,r_sat_ecef,...
        R_eci2ecef,w_eci);
    
    % partial derivatives of pseudorange rate for deputy
    [drhodotdr_d,drhodotdv_d] = dydX_pseudorange_rate(rd_eci,r_sat_eci,...
        vd_eci,v_sat_eci);
    
    % ------------------------------
    % Assemble measurement Jacobian.
    % ------------------------------
    
    % preallocates matrix for measurement Jacobian
    C = zeros(4*l,16);
    
    % storing partial derivatives of chief pseudoranges
    C(1:l,1:3) = drhodr_c;
    C(1:l,4:6) = drhodv_c;
    
    % storing partial derivatives of chief pseudorange rates
    C((l+1):(2*l),1:3) = drhodotdr_c;
    C((l+1):(2*l),4:6) = drhodotdv_c;
    
    % storing partial derivatives of deputy pseudoranges
    C((2*l+1):(3*l),9:11) = drhodr_d;
    C((2*l+1):(3*l),12:14) = drhodv_d;
    
    % storing partial derivatives of deputy pseudorange rates
    C((3*l+1):(4*l),9:11) = drhodotdr_d;
    C((3*l+1):(4*l),12:14) = drhodotdv_d;
    
    % storing partial derivatives of chief pseudoranges w.r.t. chief clock
    % bias
    C(1:l,7) = 1;
    
    % storing partial derivatives of chief pseudorange rates w.r.t. chief
    % clock bias drift rates
    C((l+1):(2*l),8) = 1;

    % storing partial derivatives of deputy pseudoranges w.r.t. deputy
    % clock bias
    C((2*l+1):(3*l),15) = 1;
    
    % storing partial derivatives of deputy pseudorange rates w.r.t. deputy
    % clock bias drift rates
    C((3*l+1):(4*l),16) = 1;

end