%--------------------------------------------------------------------------
%
% GPS_Constellation  MATLAB class defining a GPS constellation.
%
% Author: Josh Geiser and Tamas Kis
% Last Update: 2022-03-14
%
%--------------------------------------------------------------------------

classdef GPS_Constellation < handle
    
    properties
        gps_ephem_filepath
        ephem
        MJD_GPS         % (K×1 double) GPS time [MJD]
        GPS_week        % (K×1 double) GPS seconds of week [s]
        GPS_second      % (K×1 double) GPS seconds of week [m]
        ECEF_position   % (M×K×3 double) ECEF positions of M satellites at K time steps [m]
        ECEF_velocity   % (M×K×3 double) ECEF positions of M satellites at K time steps [m/s]
        ECI_position    % (M×K×3 double) ECI positions of M satellites at K time steps [m]
        ECI_velocity    % (M×K×3 double) ECI velocities of M satellites at K time steps [m/s]
        clock_bias      % (M×K double) clock biases of M satellites at K time steps
        clock_drift     % (M×K double) clock bias drift rates of M satellites at K time steps
        int_ambiguity   % (M×2 double) integer ambiguities of M satellites between deputy and chief
    end
    
    
    
    

    methods
        
        
        
        function obj = GPS_Constellation(simdata,n)
            %==============================================================
            % obj = GPS_Constellation(simdata)
            %
            % Constructor.
            %--------------------------------------------------------------
            %
            % ------
            % INPUT:
            % ------
            %   simdata - (1×1 struct) simulation data from an orbit
            %             simulation
            %   n       - (1×1 double) number of spacecraft in swarm
            %
            % -------
            % OUTPUT:
            % -------
            %   obj     - (1×1 GPS_Constellation) GPS constellation object
            %
            %==============================================================
            
            obj.gps_ephem_filepath = 'ephem.csv';
            obj.ephem = obj.get_CSV_data(obj.gps_ephem_filepath);
            obj.MJD_GPS = simdata.MJD_GPS;
            [obj.GPS_week,obj.GPS_second] = gps2wks(obj.MJD_GPS);
            
            % defaults "n" to 1 if not input
            if (nargin < 2) || isempty(n)
                n = 1;
            end

            % initialize GPS satellite states
            obj.init_satellite_states(simdata,n);

        end
        
        
        
        function rho = get_pseudorange(obj,r_rcv_ecef,b_rcv,SVID,...
                MJD_GPS,iono_delay,noise)
            %==============================================================
            % rho = get_pseudorange(obj,r_rcv_ecef,b_rcv,SVID,MJD_GPS)
            % rho = get_pseudorange(obj,r_rcv_ecef,b_rcv,SVID,MJD_GPS,...
            %   iono_delay,noise)
            %
            % Pseudorange measurement(s) from specified GPS satellite(s).
            %--------------------------------------------------------------
            %
            % ------
            % INPUT:
            % ------
            %   obj         - (1×1 GPS_Constellation) object of this class
            %   r_rcv_ecef  - (3×1 double) receiver position resolved in
            %                 ECEF frame [m]
            %   b_rcv       - (1×1 double) receiver clock bias [m]
            %   SVID        - (n×1 double) SVIDs of "n" GPS satellites [-]
            %   MJD_GPS     - (1×1 double) current GPS time [MJD]
            %   iono_delay  - (1×1 logical) (OPTIONAL) "true" if 
            %                 ionospheric path delay should be included, 
            %                 "false" otherwise (defaults to "false")
            %   noise       - (1×1 logical) (OPTIONAL) "true" if random 
            %                 noise should be included, "false" otherwise 
            %                 (defaults to false)
            %
            % -------
            % OUTPUT:
            % -------
            %   rho         - (n×1 double) pseudorange measurement(s) from
            %                 specified GPS satellite(s) [m]
            %
            %==============================================================

            % defaults "iono_delay" to "false" if not input
            if (nargin < 6) || isempty(iono_delay)
                iono_delay = false;
            end
        
            % defaults "noise" to "false" if not input
            if (nargin < 7) || isempty(noise)
                noise = false;
            end

            % number of GPS satellites to get pseudoranges from
            n = length(SVID);

            % preallocates vector to store pseudorange measurements
            rho = zeros(n,1);

            % ECEF positions [m] and clock biases [m] of specified GPS 
            % satellites
            r_sat_ecef = obj.get_ECEF_position(SVID,MJD_GPS);
            b_sat = obj.get_clock_bias(SVID,MJD_GPS);

            % pseudorange measurements [m]
            for k = 1:n
                rho(k) = pseudorange(r_rcv_ecef,r_sat_ecef(:,k),b_rcv,...
                    b_sat(k),MJD_GPS,iono_delay,noise);
            end

        end



        function rho_dot = get_pseudorange_rate(obj,r_rcv_ecef,...
                v_rcv_ecef,b_dot_rcv,SVID,MJD_GPS,noise)
            %==============================================================
            % rho_dot = get_pseudorange_rate(obj,r_rcv_ecef,v_rcv_ecef,...
            %   b_dot_rcv,SVID,MJD_GPS)
            % rho_dot = get_pseudorange_rate(obj,r_rcv_ecef,v_rcv_ecef,...
            %   b_dot_rcv,SVID,MJD_GPS,noise)
            %
            % Pseudorange rate measurement(s) from specified GPS 
            % satellite(s).
            %--------------------------------------------------------------
            %
            % ------
            % INPUT:
            % ------
            %   obj         - (1×1 GPS_Constellation) object of this class
            %   r_rcv_ecef  - (3×1 double) receiver position resolved in
            %                 ECEF frame [m]
            %   v_rcv_ecef  - (3×1 double) receiver ECEF velocity resolved
            %                 in ECEF frame [m/s]
            %   b_dot_rcv   - (1×1 double) receiver clock bias drift rate
            %                 [m/s]
            %   SVID        - (n×1 double) SVIDs of "n" GPS satellites [-]
            %   MJD_GPS     - (1×1 double) current GPS time [MJD]
            %   noise       - (1×1 logical) (OPTIONAL) "true" if random 
            %                 noise should be included, "false" otherwise 
            %                 (defaults to false)
            %
            % -------
            % OUTPUT:
            % -------
            %   rho_dot     - (n×1 double) pseudorange rate measurement(s) 
            %                 from specified GPS satellite(s) [m]
            %
            %==============================================================
        
            % defaults "noise" to "false" if not input
            if (nargin < 7) || isempty(noise)
                noise = false;
            end

            % number of GPS satellites to get pseudorange rates from
            n = length(SVID);

            % preallocates vector to store pseudorange rate measurements
            rho_dot = zeros(n,1);

            % ECEF positions [m], ECEF velocities [m/s], and clock bias
            % drift rates [m/s] of specified GPS satellites
            r_sat_ecef = obj.get_ECEF_position(SVID,MJD_GPS);
            v_sat_ecef = obj.get_ECEF_velocity(SVID,MJD_GPS);
            b_dot_sat = obj.get_clock_drift(SVID,MJD_GPS);

            % pseudorange rate measurements [m/s]
            for k = 1:n
                rho_dot(k) = pseudorange_rate(r_rcv_ecef,...
                    r_sat_ecef(:,k),v_rcv_ecef,v_sat_ecef(:,k),...
                    b_dot_rcv,b_dot_sat(k),noise);
            end

        end



        function rho_phi = get_carrier_phase(obj,r_rcv_ecef,b_rcv,SVID,...
                MJD_GPS,N,iono_delay,noise,RCVID)
            %==============================================================
            % rho_phi = get_carrier_phase(obj,r_rcv_ecef,b_rcv,SVID,...
            %   MJD_GPS,N)
            % rho_phi = get_carrier_phase(obj,r_rcv_ecef,b_rcv,SVID,...
            %   MJD_GPS,N,iono_delay,noise,RCVID)
            %
            % Carrier phase measurement(s) from specified GPS satellite(s).
            %--------------------------------------------------------------
            %
            % ------
            % INPUT:
            % ------
            %   obj         - (1×1 GPS_Constellation) object of this class
            %   r_rcv_ecef  - (3×1 double) receiver position resolved in
            %                 ECEF frame [m]
            %   b_rcv       - (1×1 double) receiver clock bias [m]
            %   SVID        - (n×1 double) SVIDs of "n" GPS satellites [-]
            %   MJD_GPS     - (1×1 double) current GPS time [MJD]
            %   N           - (n×1 double) (OPTIONAL) integer ambiguities
            %                 of the specified GPS satellites [-]
            %   iono_delay  - (1×1 logical) (OPTIONAL) "true" if 
            %                 ionospheric path delay should be included, 
            %                 "false" otherwise (defaults to "false")
            %   noise       - (1×1 logical) (OPTIONAL) "true" if random 
            %                 noise should be included, "false" otherwise 
            %                 (defaults to false)
            %   RCVID       - (1×1 double) receiver ID [-]
            %
            % -------
            % OUTPUT:
            % -------
            %   rho_phi     - (n×1 double) carrier phase measurement(s) 
            %                 from specified GPS satellite(s) [m]
            %
            %==============================================================
            
            % # of GPS satellites to get carrier phase measurements from
            l = length(SVID);

            % defaults receiver ID to 1 if not input
            if (nargin < 9) || isempty(RCVID)
                RCVID = 1;
            end
            
            % defaults "N" to the true integer ambiguities if not input
            if (nargin < 6) || isempty(N)
                N = zeros(l,1);
                for k = 1:l
                    N(k) = obj.int_ambiguity(k,RCVID);
                end
            end
            
            % defaults "iono_delay" to "false" if not input
            if (nargin < 7) || isempty(iono_delay)
                iono_delay = false;
            end
            
            % defaults "noise" to "false" if not input
            if (nargin < 8) || isempty(noise)
                noise = false;
            end
            
            % preallocates vector to store carrier phase measurements
            rho_phi = zeros(l,1);
            
            % ECEF positions [m] and clock biases [m] of specified GPS 
            % satellites
            r_sat_ecef = obj.get_ECEF_position(SVID,MJD_GPS);
            b_sat = obj.get_clock_bias(SVID,MJD_GPS);
            
            % carrier phase measurements [m]
            for k = 1:l
                rho_phi(k) = carrier_phase(r_rcv_ecef,r_sat_ecef(:,k),...
                    b_rcv,b_sat(k),MJD_GPS,N(k),iono_delay,noise);
            end
            
        end
        


        function [rho_gr,rho,rho_phi] = get_graphic(obj,r_rcv_ecef,...
                b_rcv,SVID,MJD_GPS,N,iono_delay,noise,RCVID)
            %==============================================================
            % rho_gr = get_graphic(obj,r_rcv_ecef,b_rcv,SVID,MJD_GPS,N)
            % rho_gr = get_graphic(obj,r_rcv_ecef,b_rcv,SVID,MJD_GPS,N,...
            %   iono_delay,noise,RCVID)
            % [rho_gr,rho,rho_phi] = get_graphic(__)
            %
            % Group and phase ionospheric correction (GRAPHIC)
            % measurement(s) from specified GPS satellite(s).
            %--------------------------------------------------------------
            %
            % ------
            % INPUT:
            % ------
            %   obj         - (1×1 GPS_Constellation) object of this class
            %   r_rcv_ecef  - (3×1 double) receiver position resolved in
            %                 ECEF frame [m]
            %   b_rcv       - (1×1 double) receiver clock bias [m]
            %   SVID        - (n×1 double) SVIDs of "n" GPS satellites [-]
            %   MJD_GPS     - (1×1 double) current GPS time [MJD]
            %   N           - (n×1 double) (OPTIONAL) integer ambiguities
            %                 of the specified GPS satellites [-]
            %   iono_delay  - (1×1 logical) (OPTIONAL) "true" if 
            %                 ionospheric path delay should be included, 
            %                 "false" otherwise (defaults to "false")
            %   noise       - (1×1 logical) (OPTIONAL) "true" if random 
            %                 noise should be included, "false" otherwise 
            %                 (defaults to false)
            %   RCVID       - (1×1 double) (OPTIONAL) receiver ID [-]
            %
            % -------
            % OUTPUT:
            % -------
            %   rho_gr      - (n×1 double) GRAPHIC measurement(s) from
            %                 specified GPS satellite(s) [m]
            %   rho         - (n×1 double) pseudorange measurements from
            %                 specified GPS satellite(s) [m]
            %   rho_phi     - (n×1 double) carrier phase measurements from
            %                 specified GPS satellite(s) [m]
            %
            %==============================================================
            
            % number of GPS satellites to get GRAPHIC measurements from
            n = length(SVID);
                        
            % defaults receiver ID to 1 if not input
            if (nargin < 9) || isempty(RCVID)
                RCVID = 1;
            end

            % defaults "N" to the true integer ambiguities if not input
            if (nargin < 6) || isempty(N)
                N = zeros(n,1);
                for k = 1:n
                    N(k) = obj.int_ambiguity(k,RCVID);
                end
            end
            
            % defaults "iono_delay" to "false" if not input
            if (nargin < 7) || isempty(iono_delay)
                iono_delay = false;
            end
            
            % defaults "noise" to "false" if not input
            if (nargin < 8) || isempty(noise)
                noise = false;
            end
            
            % preallocates vectors to store GRAPHIC, pseudorange, and
            % carrier phase measurements
            rho_gr = zeros(n,1);
            rho = zeros(n,1);
            rho_phi = zeros(n,1);
            
            % ECEF positions [m] and clock biases [m] of specified GPS 
            % satellites
            r_sat_ecef = obj.get_ECEF_position(SVID,MJD_GPS);
            b_sat = obj.get_clock_bias(SVID,MJD_GPS);
            
            % GRAPHIC, pseudorange, and carrier phase measurements [m]
            for k = 1:n
                [rho_gr(k),rho(k),rho_phi(k)] = graphic(r_rcv_ecef,...
                    r_sat_ecef(:,k),b_rcv,b_sat(k),MJD_GPS,N(k),...
                    iono_delay,noise);
            end
            
        end
        
        

        function [rho_sdcp,rho_phi_d,rho_phi_c] = get_sdcp(obj,rd_ecef,...
                rc_ecef,bd,bc,SVID,MJD_GPS,Nd,Nc,iono_delay,noise,...
                RCVID_d,RCVID_c)
            %==============================================================
            % rho_sdcp = get_sdcp(obj,rd_ecef,rc_ecef,bd,bc,SVID,MJD_GPS)
            % rho_sdcp = get_sdcp(obj,rd_ecef,rc_ecef,bd,bc,SVID,...
            %   MJD_GPS,Nd,Nc,iono_delay,noise,RCVID_d,RCVID_c)
            % [rho_sdcp,rho_phi_d,rho_phi_c] = get_sdcp(__)
            %
            % Single-difference carrier phase (SDCP) measurement(s) from 
            % specified GPS satellite(s).
            %--------------------------------------------------------------
            %
            % ------
            % INPUT:
            % ------
            %   obj         - (1×1 GPS_Constellation) object of this class
            %   rd_ecef     - (3×1 double) deputy position resolved in ECEF
            %                 frame [m]
            %   rc_ecef     - (3×1 double) chief position resolved in ECEF
            %                 frame [m]
            %   bd          - (1×1 double) deputy clock bias [m]
            %   bc          - (1×1 double) chief clock bias [m]
            %   SVID        - (n×1 double) SVIDs of "n" GPS satellites [-]
            %   MJD_GPS     - (1×1 double) current GPS time [MJD]
            %   Nd          - (n×1 double) (OPTIONAL) integer ambiguities 
            %                 between deputy and specified GPS satellite(s)
            %                 [-]
            %   Nc          - (n×1 double) (OPTIONAL) integer ambiguities 
            %                 between chief and specified GPS satellite(s)
            %                 [-]
            %   iono_delay  - (1×1 logical) (OPTIONAL) "true" if 
            %                 ionospheric path delay should be included, 
            %                 "false" otherwise (defaults to "false")
            %   noise       - (1×1 logical) (OPTIONAL) "true" if random 
            %                 noise should be included, "false" otherwise 
            %                 (defaults to false)
            %   RCVID_d     - (1×1 double) deputy receiver ID [-]
            %   RCVID_c     - (1×1 double) chief receiver ID [-]
            %
            % -------
            % OUTPUT:
            % -------
            %   rho_sdcp    - (n×1 double) SDCP measurement(s) from
            %                 specified GPS satellite(s) [m]
            %   rho_phi_d   - (n×1 double) deputy carrier phase 
            %                 measurement(s) from specified GPS 
            %                 satellite(s) [m]
            %   rho_phi_c   - (n×1 double) chief carrier phase 
            %                 measurement(s) from specified GPS 
            %                 satellite(s) [m]
            %
            %==============================================================

            % number of GPS satellites to get SDCP measurements from
            l = length(SVID);
            
            % defaults deputy receiver ID to 2 if not input
            if (nargin < 12) || isempty(RCVID_d)
                RCVID_d = 2;
            end

            % defaults chief receiver ID to 1 if not input
            if (nargin < 13) || isempty(RCVID_c)
                RCVID_c = 1;
            end

            % defaults "Nd" to the true integer ambiguities between GPS
            % satellites and deputy if not input
            if (nargin < 8) || isempty(Nd)
                Nd = zeros(l,1);
                for k = 1:l
                    Nd(k) = obj.int_ambiguity(k,RCVID_d);
                end
            end

            % defaults "Nc" to the true integer ambiguities between GPS
            % satellites and chief if not input
            if (nargin < 9) || isempty(Nc)
                Nc = zeros(l,1);
                for k = 1:l
                    Nc(k) = obj.int_ambiguity(k,RCVID_c);
                end
            end
            
            % defaults "iono_delay" to "false" if not input
            if (nargin < 10) || isempty(iono_delay)
                iono_delay = false;
            end
            
            % defaults "noise" to "false" if not input
            if (nargin < 11) || isempty(noise)
                noise = false;
            end
            
            % preallocates vectors to store SDCP and carrier phase 
            % measurements
            rho_sdcp = zeros(l,1);
            rho_phi_d = zeros(l,1);
            rho_phi_c = zeros(l,1);
            
            % ECEF positions [m] and clock biases [m] of specified GPS 
            % satellites
            r_sat_ecef = obj.get_ECEF_position(SVID,MJD_GPS);
            b_sat = obj.get_clock_bias(SVID,MJD_GPS);
            
            % SDCP measurements [m]
            for k = 1:l
                [rho_sdcp(k),rho_phi_d(k),rho_phi_c(k)] = sdcp(rd_ecef,...
                    rc_ecef,r_sat_ecef(:,k),bd,bc,b_sat(k),MJD_GPS,...
                    Nd(k),Nc(k),iono_delay,noise);
            end
            
        end
        
        
        
        function SVID = get_closest_SVIDs(obj,r_rcv_ecef,MJD_GPS,l)
            %==============================================================
            % SVID = get_closest_SVIDs(obj,r_rcv_ecef,MJD_GPS,l)
            %
            % SVIDs of the "l" closest GPS satellites to a receiver.
            %--------------------------------------------------------------
            %
            % ------
            % INPUT:
            % ------
            %   obj         - (1×1 GPS_Constellation) object of this class
            %   r_rcv_ecef  - (3×1 double) receiver position resolved in
            %                 ECEF frame [m]
            %   MJD_GPS     - (1×1 double) current GPS time [MJD]
            %   l           - (1×1 double) (OPTIONAL) number of SVIDs to
            %                 return (defaults to 4)
            %
            % -------
            % OUTPUT:
            % -------
            %   SVID        - (n×1 double) SVIDs of the "n" closest GPS
            %                 satellites
            %
            %==============================================================
            
            % total number of GPS satellites
            M = max(size(obj.ephem));
            
            % preallocates vector to true store ranges
            ranges = zeros(M,1);
            
            % index corresponding to current sample time
            i = find(MJD_GPS == obj.MJD_GPS);
            
            % calculates ranges to all GPS satellites
            for k = 1:M
                r_sat_ecef = squeeze(obj.ECEF_position(k,i,1:3));
                ranges(k) = inorm(r_sat_ecef-r_rcv_ecef);
            end
            
            % SVIDs of closest "n" GPS satellites
            [~,SVID] = mink(ranges,l);

            % reshapes vector storing SVIDs
            SVID = SVID.';

        end
        
        
        
        function r_sat_ecef = get_ECEF_position(obj,SVID,MJD_GPS)
            %==============================================================
            % r_sat_ecef = get_ECEF_position(obj,SVID,MJD_GPS)
            %
            % GPS satellite ECEF position.
            %--------------------------------------------------------------
            %
            % ------
            % INPUT:
            % ------
            %   obj         - (1×1 GPS_Constellation) object of this class
            %   SVID        - (n×1 double) SVIDs of "n" GPS satellites [-]
            %   MJD_GPS     - (1×1 double) GPS time [MJD]
            %
            % -------
            % OUTPUT:
            % -------
            %   r_sat_ecef  - (3×n double) ECEF position(s) of specified
            %                 GPS satellite(s) at specified time [m]
            %
            %==============================================================
            
            % index corresponding to current GPS time
            i = find(MJD_GPS == obj.MJD_GPS);
            
            % number of positions to return
            n = length(SVID);
            
            % preallocates array to store positions
            r_sat_ecef = zeros(3,n);
            
            % ECEF position(s) of specified GPS satellite(s) [m]
            for k = 1:n
                r_sat_ecef(:,k) = squeeze(obj.ECEF_position(SVID(k),i,:));
            end
            
        end
        
        
        
        function v_sat_ecef = get_ECEF_velocity(obj,SVID,MJD_GPS)
            %==============================================================
            % v_sat_ecef = get_ECEF_velocity(obj,SVID,MJD_GPS)
            %
            % GPS satellite ECEF velocity.
            %--------------------------------------------------------------
            %
            % ------
            % INPUT:
            % ------
            %   obj         - (1×1 GPS_Constellation) object of this class
            %   SVID        - (n×1 double) SVIDs of "n" GPS satellites [-]
            %   MJD_GPS     - (1×1 double) GPS time [MJD]
            %
            % -------
            % OUTPUT:
            % -------
            %   v_sat_ecef  - (3×n double) ECEF velocity(ies) of specified 
            %                 GPS satellite(s) at specified time [m/s]
            %
            %==============================================================
            
            % index corresponding to current GPS time
            i = find(MJD_GPS == obj.MJD_GPS);
            
            % number of velocities to return
            n = length(SVID);

            % preallocates array to store velocities
            v_sat_ecef = zeros(3,n);

            % ECEF velocity(ies) of specified GPS satellite(s) [m/s]
            for k = 1:n
                v_sat_ecef(:,k) = squeeze(obj.ECEF_velocity(SVID(k),i,:));
            end
            
        end



        function r_sat_eci = get_ECI_position(obj,SVID,MJD_GPS)
            %==============================================================
            % r_sat_eci = get_ECI_position(obj,SVID,MJD_GPS)
            %
            % GPS satellite ECI position.
            %--------------------------------------------------------------
            %
            % ------
            % INPUT:
            % ------
            %   obj         - (1×1 GPS_Constellation) object of this class
            %   SVID        - (n×1 double) SVIDs of "n" GPS satellites [-]
            %   MJD_GPS     - (1×1 double) GPS time [MJD]
            %
            % -------
            % OUTPUT:
            % -------
            %   r_sat_eci   - (3×n double) ECI position(s) of specified GPS
            %                 satellite(s) at specified time [m]
            %
            %==============================================================
            
            % index corresponding to current GPS time
            i = find(MJD_GPS == obj.MJD_GPS);

            % number of positions to return
            n = length(SVID);
            
            % preallocates array to store positions
            r_sat_eci = zeros(3,n);
            
            % ECI position(s) of specified GPS satellite(s) [m]
            for k = 1:n
                r_sat_eci(:,k) = squeeze(obj.ECI_position(SVID(k),i,:));
            end
            
        end
        
        

        function v_sat_eci = get_ECI_velocity(obj,SVID,MJD_GPS)
            %==============================================================
            % v_sat_eci = get_ECI_velocity(obj,SVID,MJD_GPS)
            %
            % GPS satellite ECI velocity.
            %--------------------------------------------------------------
            %
            % ------
            % INPUT:
            % ------
            %   obj         - (1×1 GPS_Constellation) object of this class
            %   SVID        - (n×1 double) SVIDs of "n" GPS satellites [-]
            %   MJD_GPS     - (1×1 double) GPS time [MJD]
            %
            % -------
            % OUTPUT:
            % -------
            %   v_sat_eci   - (3×n double) ECI velocity(ies) of specified 
            %                 GPS satellite(s) at specified time [m/s]
            %
            %==============================================================
            
            % index corresponding to current GPS time
            i = find(MJD_GPS == obj.MJD_GPS);
            
            % number of velocities to return
            n = length(SVID);

            % preallocates array to store velocities
            v_sat_eci = zeros(3,n);

            % ECI velocity(ies) of specified GPS satellite(s) [m/s]
            for k = 1:n
                v_sat_eci(:,k) = squeeze(obj.ECI_velocity(SVID(k),i,:));
            end
            
        end



        function b_sat = get_clock_bias(obj,SVID,MJD_GPS)
            %==============================================================
            % b_sat = get_clock_bias(obj,SVID,MJD_GPS)
            %
            % GPS satellite clock bias.
            %--------------------------------------------------------------
            %
            % ------
            % INPUT:
            % ------
            %   obj         - (1×1 GPS_Constellation) object of this class
            %   SVID    - (n×1 double) SVIDs of "n" GPS satellites [-]
            %   MJD_GPS - (1×1 double) GPS time [MJD]
            %
            % -------
            % OUTPUT:
            % -------
            %   b_sat   - (n×1 double) clock bias(es) of specified GPS 
            %             satellite(s) at specified time [m]
            %
            %==============================================================
            
            % index corresponding to current GPS time
            i = find(MJD_GPS == obj.MJD_GPS);
            
            % number of clock biases to return
            n = length(SVID);

            % preallocates array to store clock biases
            b_sat = zeros(n,1);
            
            % clock bias(es) of specified GPS satellite(s) [m]
            for k = 1:n
                b_sat(k) = obj.clock_bias(SVID(k),i);
            end
            
        end
        
        
        
        function b_dot_sat = get_clock_drift(obj,SVID,MJD_GPS)
            %==============================================================
            % b_dot_sat = get_clock_drift(obj,SVID,MJD_GPS)
            %
            % GPS satellite clock bias drift rate.
            %--------------------------------------------------------------
            %
            % ------
            % INPUT:
            % ------
            %   obj         - (1×1 GPS_Constellation) object of this class
            %   SVID        - (n×1 double) SVIDs of "n" GPS satellites [-]
            %   MJD_GPS     - (1×1 double) GPS time [MJD]
            %
            % -------
            % OUTPUT:
            % -------
            %   b_dot_sat   - (1×1 double) clock bias drift rate of 
            %                 specified GPS satellite at specified time 
            %                 [m/s]
            %
            %==============================================================
            
            % index corresponding to current GPS time
            i = find(MJD_GPS == obj.MJD_GPS);
            
            % number of clock bias drift rates to return
            n = length(SVID);
            
            % preallocates array to store clock bias drift rates
            b_dot_sat = zeros(n,1);

            % clock bias drift rate(s) of specified GPS satellite(s) [m/s]
            for k = 1:n
                b_dot_sat(k) = obj.clock_drift(SVID(k),i);
            end
            
        end


        
        function v_sat_ecef = compute_ECEF_velocity(obj,SVID,GPS_s)
            %==============================================================
            % r_sat_ecef = compute_ECEF_velocity(obj,SVID,GPS_s)
            %
            % ECEF position of the GPS satellite with the specified SVID at 
            % the specified time.
            %--------------------------------------------------------------
            %
            % ------
            % INPUT:
            % ------
            %   obj         - (1×1 GPS_Constellation) object of this class
            %   SVID        - (1×1 double) GPS satellite vehicle ID
            %   GPS_s       - (1×1 double) GPS seconds of GPS week [s]
            %
            % -------
            % OUTPUT:
            % -------
            %   v_sat_ecef  - (3×1 double) ECEF velocity of specified GPS 
            %                 satellite [m/s]
            %
            %==============================================================
            
            % time step to use for forward difference approximation [s]
            dt = 0.001;

            % ECEF position of specified GPS satellite at current time [m]
            r_sat_ecef = obj.ephem2ecef(obj.ephem(SVID,:),GPS_s);

            % ECEF position of specified GPS satellite "dt" seconds in the
            % future [m]
            r_sat_ecef_dt = obj.ephem2ecef(obj.ephem(SVID,:),GPS_s+dt);

            % (approximate) ECEF velocity of specified GPS satellite
            % resolved in ECEF frame [m/s]
            v_sat_ecef = (r_sat_ecef_dt-r_sat_ecef)/dt;
            
        end



        function init_satellite_states(obj,simdata,n)
            %==============================================================
            % init_satellite_state(obj,simdata,n)
            %
            % Initializes "ECEF_position" and "clock_bias" fields.
            %==============================================================
            
            % number of time steps
            K = length(obj.GPS_second);
            
            % number of GPS satellites
            M = max(size(obj.ephem));
            
            % preallocates M×K×3 array to store ECEF positions
            obj.ECEF_position = zeros(M,K,3);
            
            % preallocates M×K array to store clock biases
            obj.clock_bias = zeros(M,K,1);
            
            % GPS constellation data at every sample time
            for i = 1:M
                for j = 1:K
                    
                    % computes ECEF position and clock bias from ephemeris 
                    % data [m]
                    [r_sat_ecef,b_sat] = obj.ephem2ecef(obj.ephem(i,:),...
                        obj.GPS_second(j));
                    
                    % computes ECEF velocity [m/s]
                    v_sat_ecef = compute_ECEF_velocity(obj,i,...
                        obj.GPS_second(j));
                    
                    % rotation matrix (ECEF --> ECI)
                    R_ecef2eci = simdata.R_ecef2eci(:,:,j);
                    w_eci = simdata.w_eci(:,j);
                    
                    % ECI position [m] and velocity [m/s]
                    [r_sat_eci,v_sat_eci] = ecef2eci(r_sat_ecef,...
                        v_sat_ecef,w_eci,R_ecef2eci);
                    
                    % clock bias drift rate [m/s]
                    b_dot_sat = obj.ephem(i,:).SVclockDrift;
                    
                    % stores data
                    obj.ECEF_position(i,j,:) = r_sat_ecef;
                    obj.ECEF_velocity(i,j,:) = v_sat_ecef;
                    obj.ECI_position(i,j,:) = r_sat_eci;
                    obj.ECI_velocity(i,j,:) = v_sat_eci;
                    obj.clock_bias(i,j) = b_sat;
                    obj.clock_drift(i,j) = b_dot_sat;
                    
                end
            end

            % randomly picks integer ambiguities between 50 and 200 for 
            % each spacecraft in a swarm with each GPS satellite
            obj.int_ambiguity = rand2(50,200,[M,n],'int');
            
        end
        
        
        
        function [ECEF_est,delta_t_SV] = ephem2ecef(~,ephem,tx_time)
            %==============================================================
            % [ECEF_est,delta_t_SV] = ephem2ecef(~,ephem,tx_time)
            %
            % ECEF position and clock bias from ephemeris data and
            % transmission time.
            %--------------------------------------------------------------
            %
            % ------
            % INPUT:
            % ------
            %   ephem   - (1×1 struct) GPS ephemeris data
            %   tx_time - (1×1 double) transmission time [s]
            %
            %==============================================================

            % Define some constants
            MU = 3.986005e14;
            OM_E_DOT = 7.2921151467e-5;
            C = 2.99792458e8;
            F = -4.442807633e-10;

            % 19-step procedure for calculating satellite ECEF position from ephemeris
            a = ephem.sqrtA ^ 2;
            n = sqrt(MU/(a^3)) + ephem.DeltaN;
            t_k = tx_time - ephem.Toe;
            M_k = ephem.M0 + n*t_k;
            E_k = M2E(M_k, ephem.Eccentricity);
            sin_v_k = sqrt(1 - ephem.Eccentricity^2)*sin(E_k) / (1 - ephem.Eccentricity*cos(E_k));
            cos_v_k = (cos(E_k)-ephem.Eccentricity) / (1 - ephem.Eccentricity*cos(E_k)); 
            v_k = atan2(sin_v_k, cos_v_k);
            phi_k = v_k + ephem.omega;
            delta_phi_k = ephem.Cus * sin(2*phi_k) + ephem.Cuc * cos(2*phi_k);
            u_k = phi_k + delta_phi_k;
            delta_r_k = ephem.Crs * sin(2*phi_k) + ephem.Crc * cos(2*phi_k);
            delta_i_k = ephem.Cis * sin(2*phi_k) + ephem.Cic * cos(2*phi_k); 
            Om_k = ephem.Omega0 - OM_E_DOT*tx_time + ephem.OmegaDot*t_k; % + -0.00000575;
            r_k = a * (1 - ephem.Eccentricity * cos(E_k)) + delta_r_k;
            i_k = ephem.Io + ephem.IDOT*t_k + delta_i_k;
            x_p = r_k * cos(u_k);
            y_p = r_k * sin(u_k); 
            X_ECEF = x_p * cos(Om_k) - y_p * cos(i_k) * sin(Om_k);
            Y_ECEF = x_p * sin(Om_k) + y_p * cos(i_k) * cos(Om_k);
            Z_ECEF = y_p * sin(i_k);
            ECEF_est = [X_ECEF; Y_ECEF; Z_ECEF];

            % Calculation of satellite's clock bias term
            delta_t_r = F * ephem.Eccentricity * ephem.sqrtA * sin(E_k);
            delta_t_SV = ephem.SVclockBias ...
                         + ephem.SVclockDrift * (tx_time - ephem.TransTime) ...
                         + ephem.SVclockDriftRate * (tx_time - ephem.TransTime)^2 ...
                         + delta_t_r ...
                         - ephem.TGD;
            delta_t_SV = delta_t_SV * C;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Helper function for reading CSV file into nice data table format
        function data = get_CSV_data(~,filename)
            data = readtable(filename);
        end
        
    end 
end