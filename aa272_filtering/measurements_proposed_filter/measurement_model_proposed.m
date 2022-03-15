%==========================================================================
%
% measurement_model_proposed  Measurement model. 
%
%   y = measurement_model_proposed(x,t,iono_delay,noise,prop,gps)
%   y = measurement_model_proposed(x,t,iono_delay,noise,prop,gps,l)
%   [y,rho_c,rho_d,SVID] = measurement_model_proposed(__)
%
% Author: Tamas Kis
% Last Update: 2022-03-14
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x           - (16×1 double) state vector
%                   -->   1-3. rc_eci - chief position resolved in ECI
%                                       frame [m]
%                   -->   4-6. vc_eci - chief inertial velocity resolved in
%                                       ECI frame [m/s]
%                   -->     7. bc     - chief clock bias [m]
%                   -->     8. bc_dot - chief clock bias drift rate [m/s]
%                   -->  9-11. rd_eci - deputy position resolved in ECI
%                                       frame [m]
%                   --> 11-14. vd_eci - deputy inertial velocity resolved
%                                       in ECI frame [m/s]
%                   -->    15. bd     - deputy clock bias [m]
%                   -->    16. bd_dot - deputy clock bias drift rate [m/s]
%   t           - (1×1 double) simulation time [s]
%   iono_delay  - (1×1 logical) "true" if ionospheric path delay should be 
%                 included, "false" otherwise
%   noise       - (1×1 logical) "true" if true measurement noise should be 
%                 included, "false" otherwise
%   prop        - (1×1 struct) propagator settings (see 
%                 "initialize_propagator" function for full definition)
%   gps         - (1×1 GPS_Constellation) GPS constellation object
%   l           - (1×1 double) (OPTIONAL) number of GNSS satellites used
%                 for measurements (defaults to 4)
%   Nc          - (l×1 double) (OPTIONAL) integer ambiguities between
%                 chief and "l" closest satellites to chief
%   Nd          - (l×1 double) (OPTIONAL) integer ambiguities between
%                 deputy and "l" closest satellites to chief
%
% -------
% OUTPUT:
% -------
%   y       - (4l×1 double) measurement
%               -->       1-l. "l" GRAPHIC Hatch measurements for chief 
%                              satellite
%               -->  (l+1)-2l. "l" pseudorange rates for chief satellite
%               --> (2l+1)-3l. "l" SDCP Hatch measurements between chief
%                              and deputy satellites
%               --> (3l+1)-4l. "l" pseudorange rates for deputy satellite
%   rho_c   - (l×1 double) chief pseudorange measurements
%   rho_d   - (l×1 double) deputy pseudorange measurements
%   SVID    - (l×1 double) SVIDs of "l" closest satellites to chief
%
%==========================================================================
function [y,rho_c,rho_d,SVID] = measurement_model_proposed(x,t,...
    iono_delay,noise,prop,gps,l,Nc,Nd)
    
    % initializes integer ambiguities to empty vector if not input
    if (nargin < 8)
        Nc = [];
        Nd = [];
    end

    % ---------------------
    % Unpacks state vector.
    % ---------------------
    
    % chief position [m] and inertial velocity [m/s] resolved in ECI frame
    rc_eci = x(1:3);
    vc_eci = x(4:6);
    
    % chief block bias [m] and clock bias drift rate [m/s]
    bc = x(7);
    bc_dot = x(8);
    
    % deputy position [m] and inertial velocity [m/s] resolved in ECI frame
    rd_eci = x(9:11);
    vd_eci = x(12:14);
    
    % deputy block bias [m] and clock bias drift rate [m/s]
    bd = x(15);
    bd_dot = x(16);
    
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
    
    % ------------
    % ECEF states.
    % ------------
    
    % chief position [m] and ECEF velocity [m/s] resolved in ECEF frame
    [rc_ecef,vc_ecef] = eci2ecef(rc_eci,vc_eci,w_eci,R_eci2ecef);
    
    % deputy position [m] and ECEF velocity [m/s] resolved in ECEF frame
    [rd_ecef,vd_ecef] = eci2ecef(rd_eci,vd_eci,w_eci,R_eci2ecef);
    
    % ----------------------------------------------
    % Pseudorange and pseudorange rate measurements.
    % ----------------------------------------------
    
    % defaults "l" to 4 if not input
    if (nargin < 7) || isempty(l)
        l = 4;
    end

    % SVIDs of the "l" closest GPS satellites to chief
    SVID = gps.get_closest_SVIDs(rc_ecef,MJD_GPS,l);
    
    % chief GRAPHIC and pseudorange measurements [m]
    [rho_gr_c,rho_c] = gps.get_graphic(rc_ecef,bc,SVID,MJD_GPS,Nc,...
        iono_delay,noise);
    
    % chief pseudorange rate measurements [m/s]
    rho_dot_c = gps.get_pseudorange_rate(rc_ecef,vc_ecef,bc_dot,SVID,...
        MJD_GPS,noise);
    
    % deputy/chief SDCP measurements [m]
    rho_sdcp_dc = gps.get_sdcp(rd_ecef,rc_ecef,bd,bc,SVID,MJD_GPS,Nd,Nc,...
        iono_delay,noise);
    
    % deputy pseudorange rate measurements [m/s]
    rho_dot_d = gps.get_pseudorange_rate(rd_ecef,vd_ecef,bd_dot,SVID,...
        MJD_GPS,noise);

    % ---------------------------------------------------------------------
    % Extra measurement (needed for Hatch filtering for true measurements).
    % ---------------------------------------------------------------------

    % deputy pseudorange measurement [m]
    rho_d = gps.get_pseudorange(rd_ecef,bd,SVID,MJD_GPS,iono_delay,noise);
    
    % ---------------------
    % Assemble measurement.
    % ---------------------
    
    y = [rho_gr_c
         rho_dot_c;
         rho_sdcp_dc;
         rho_dot_d];
    
end