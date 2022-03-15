%==========================================================================
%
% plot_filter_results  Plots the results of state estimation for a single
% state variable.
%
%   plot_filter_results(t,x)
%   plot_filter_results(t,x,x_lower,x_upper)
%   plot_filter_results(t,x,[],[],x_true)
%   plot_filter_results(t,x,x_lower,x_upper,x_true)
%   plot_filter_results(__,opts)
%
% Author: Tamas Kis
% Last Update: 2022-03-15
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   t           - (1×N double) time vector
%   x           - (1×N double) a posteriori state variable estimates
%   x_lower     - (OPTIONAL) (1×N double) lower covariance bound on state 
%                 variable estimates
%   x_upper     - (OPTIONAL) (1×N double) upper covariance bound on state 
%                 variable estimates
%   x_true      - (OPTIONAL) (1×N double) true state variable time history
%   opts        - (OPTIONAL) (struct) plot options
%       • color  - (char or 1×3 double) color scheme (defaults to 
%                  [0,0.4470,0.7410]) [rgb]
%                       --> can be specified as a name, short name, or
%                           RGB triplet
%       • digits - (1×1 double) number of digits after decimal point for
%                  scientific notation (defaults to 4)
%       • dots   - (1×1 logical) true if state estimate and true state
%                   should be plotted with dots instead of a line, false
%                   otherwise (defaults to false)
%       • error  - (1×1 logical) true if error in state estimate should be
%                  plotted, false otherwise (defaults to false)
%       • figure - (1×1 logical) true if new figure should be created,
%                   false otherwise (defaults to true)
%       • grid   - (1×1 logical) true for grid on, false otherwise
%                  (defaults to true)
%       • legend - (1×1 logical) true to include legend, false otherwise
%                  (defaults to true)
%       • M      - (1×1 double) number of standard deviations used for
%                  covariance bounds (i.e. Mσ covariance bounds) (defaults 
%                  to 1)
%       • name   - (1×1 string) name of state variable
%       • scinot - (1×1 double) 'true' or 'false' (defaults to true)
%                   --> displays numbers in scientific notation if true
%       • shaded - (1×1 logical) 'true' or 'false' (defaults to true)
%                   --> shades bounds if "true"
%                   --> draws bounds with lines if "false"
%       • symbol - (1×1 string) symbol for state variable
%       • tlabel - (1×1 logical) true if time axis label should be
%                  included, false otherwise (defaults to true)
%       • tunits - (1×1 string) time units (defaults to seconds)
%       • xlabel - (1×1 logical) true if state variable axis label should
%                  be included, false otherwise (defaults to true)
%       • xunits - (1×1 string) units of state variable
%
%==========================================================================
function plot_filter_results(t,x,x_lower,x_upper,x_true,opts)
    %TODO: DOTS VS. LINES
    %TODO: Error only
    %TODO: DOTS FOR ERROR ONLY
    %TODO: orientation of horizontal axis for fliplr for covar bounds

    % ----------------------------------------------------
    % Sets unspecified parameters to their default values.
    % ----------------------------------------------------

    % sets color scheme (defaults to the default MATLAB color)
    if (nargin < 6) || ~isfield(opts,'color')
        color = [0,0.4470,0.7410];
    else
        color = opts.color;
    end
    
    % sets number of digits after decimal point for scientific notation 
    % (defaults to 4)
    if (nargin < 6) || ~isfield(opts,'digits')
        digits = 4;
    else
        digits = opts.digits;
    end

    % sets line style (defaults to line)
    if (nargin < 6) || ~isfield(opts,'dots')
        dots = false;
    else
        dots = opts.dots;
    end

    % determines if state estimate or error in state estimate should be
    % plotted (defaults to state estimate)
    if (nargin < 6) || ~isfield(opts,'error')
        plot_error = false;
    else
        plot_error = opts.error;
    end

    % determines if new figure should be created (defaults to yes)
    if (nargin < 6) || ~isfield(opts,'figure')
        new_figure = true;
    else
        new_figure = opts.figure;
    end

    % sets grid setting (defaults to "true")
    if (nargin < 6) || ~isfield(opts,'grid')
        grid_on = true;
    else
        grid_on = opts.grid;
    end

    % sets legend setting (defaults to "true")
    if (nargin < 6) || ~isfield(opts,'legend')
        include_legend = true;
    else
        include_legend = opts.grid;
    end

    % sets # of standard deviations for covariance bounds (defaults to 1)
    if (nargin < 6) || ~isfield(opts,'M')
        M = 1;
    else
        M = opts.M;
    end

    % sets state variable name (defaults to empty string)
    if (nargin < 6) || ~isfield(opts,'name')
        name = "";
    else
        name = opts.name;
    end
    
    % sets scientific notation to be on or off (defaults to off)
    if (nargin < 6) || ~isfield(opts,'scinot')
        scinot = false;
    else
        scinot = opts.scinot;
    end

    % sets estimate bounds drawing style (defaults to shaded)
    if (nargin < 6) || ~isfield(opts,'shaded')
        shaded = true;
    else
        shaded = opts.shaded;
    end

    % sets state variable symbol (defaults to an empty string)
    if (nargin < 6) || ~isfield(opts,'symbol')
        symbol = "";
    else
        symbol = opts.symbol;
    end

    % turns time axis label on/off (defaults to on)
    if (nargin < 6) || ~isfield(opts,'tlabel')
        tlabel_on = true;
    else
        tlabel_on = opts.tlabel;
    end

    % sets time units (defaults to seconds)
    if (nargin < 6) || ~isfield(opts,'tunits')
        tunits = "s";
    else
        tunits = opts.tunits;
    end

    % turns state variable axis label on/off (defaults to on)
    if (nargin < 6) || ~isfield(opts,'xlabel')
        xlabel_on = true;
    else
        xlabel_on = opts.xlabel;
    end

    % sets state variable units (defaults to an empty string)
    if (nargin < 6) || ~isfield(opts,'xunits')
        xunits = "";
    else
        xunits = opts.xunits;
    end

    % ------
    % Plots.
    % ------

    % initializes new figure
    if new_figure
        figure('position',[540,300,700,500]);
    end

    % retain all plots
    hold on;

    % "light" color
    light_color = (1-color)*0.85+color;

    if plot_error
        x_plot = x-x_true;
        %x_plot = x;
        x_upper_plot = x_upper-x;
        x_lower_plot = x_lower-x;
    else
        x_plot = x;
        x_upper_plot = x_upper;
        x_lower_plot = x_lower;
    end

    % bounds on the state variable estimate
    if shaded
        patch([t,fliplr(t)],[x_upper_plot,fliplr(x_lower_plot)],...
            light_color,'edgecolor','none');
    else
        plot(t,x_lower_plot,'k');
        plot(t,x_upper_plot,'k','handlevisibility','off');
    end

    % plots state variable estimates
    if dots
        plot(t,x_plot,'.','markersize',7,'color',color);
    else
        plot(t,x_plot,'k','linewidth',1.5,'color',color);
    end

    % plots true state variable
    if ~plot_error
        if dots
            plot(t,x_true,'.','markersize',7,'color',color);
        else
            plot(t,x_true,'k--','linewidth',1.5,'color',color);
        end
    end

    % -------
    % Legend.
    % -------
    
    %TODO: separate function to get average bounds
    if include_legend
        
        % length of time vector
        N = length(t);
    
        % obtain average Mσ covariance bound (ignore first 10% of data due
        % to transience)
        avg_bound = mean(x_upper(round(N/10):end)-x(round(N/10):end));
    
        % converts number to string
        if scinot
            avg_bound_str = scientific_notation_string(avg_bound,digits);
        else
            avg_bound_str = ""+avg_bound;
        end
        
        % string for state estimate legend entry
        if name == ""
            estimate_str = "state estimate";
        else
            estimate_str = ""+name+" estimate";
        end
        if plot_error
            estimate_str = "error in "+estimate_str;
        end

        % string for covariance bounds legend entry
        covariance_str = "$\pm"+avg_bound_str+"\;\mathrm{"+xunits+...
                "}$ ($\pm"+M+"\sigma$)";
    
        % string for true state legend entry
        if ~plot_error
            if name ~= ""
                true_str = "true "+name;
            else
                true_str = "true state";
            end
        end

        % legend
        if plot_error
            legend(covariance_str,estimate_str,'interpreter','latex',...
                'fontsize',14,'location','northeast');
        else
            legend(covariance_str,estimate_str,true_str,'interpreter',...
                'latex','fontsize',14,'location','northeast');
        end

    end

    % ------------
    % Axis labels.
    % ------------

    % time axis label
    if tlabel_on
        xlabel("Time $\left[\mathrm{"+tunits+"}\right]$",'interpreter',...
            'latex','fontsize',18);
    end

    % state variable axis label
    if xlabel_on

        % potential name and symbol combinations
        if (name == "") && (symbol == "")
            xlabel_str = "";
        elseif (name == "") && (symbol ~= "")
            xlabel_str = symbol;
        elseif (name ~= "") && (symbol == "")
            xlabel_str = name;
        elseif (name ~= "") && (symbol ~= "")
            xlabel_str = name+", "+symbol;
        end

        % adds "error" designation to string
        if plot_error
            xlabel_str = xlabel_str+" error";
        end

        % adds units to string
        if exist('xunits','var') && (xunits ~= "")
            xlabel_str = xlabel_str+" $\left[\mathrm{"+xunits+"}\right]$";
        end

        % creates state variable axis label
        ylabel(xlabel_str,'interpreter','latex','fontsize',18);

    end
    
    % -----------------
    % Misc. formatting.
    % -----------------

    % sets time axis limits
    xlim([min(t),max(t)]);

    % sets vertical axis limits
    if plot_error
        ylim([-5*avg_bound,5*avg_bound]);
    end

    % displays grid over all other graphics objects
    if grid_on
        H = gca;
        H.Layer = 'top';
        grid on;
    end

end