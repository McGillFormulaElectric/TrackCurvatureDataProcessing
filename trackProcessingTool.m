classdef trackProcessingTool
    %TRACKPROCESSINGTOOL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % track X,Y coordinates
        X
        Y
        
        event char {mustBeMember(event,...
            {'accel', 'skidpad', 'autocross', 'endurance', ''})} = ''
        % only for endurance (optional). If no value is given assumes ~ 22km
        % track
        numberLaps
        % optional
        units char {mustBeMember(units,...
            {'m', 'ft', 'km', 'in', ''})} = ''
        
        % might have to play around with this parameter in order to get desired
        % smoothing
        % value between 0 and 1
        % lower value will give smoother curve
        butterWorthFilterParameter = 0.1 ;
        
        % track step size used for getting radii
        trackStepSize = 0.1 ;
        % track step size for output track
        trackStepSizeOut = 0.25 ;
    end
    
    methods
        function obj = set.X( obj, val )
            s = size(val) ;
            if (s(1)~=1) && (s(2)~=1)
                warning('input must be vector')
            elseif s(1) ~= 1
                obj.X = val' ;
            elseif s(2) ~= 1
                obj.X = val ;
            end
        end
        function obj = set.Y( obj, val )
            s = size(val) ;
            if (s(1)~=1) && (s(2)~=1)
                warning('input must be vector')
            elseif s(1) ~= 1
                obj.Y = val' ;
            elseif s(2) ~= 1
                obj.Y = val ;
            end
        end
        
        function track = createTrack( obj )
            % given track with x and y coordinates
            
            x = obj.X ;
            y = obj.Y ;
            
            % get distance between points
            i = 1: length(x)-1 ;
            ds = (((x(i+1)-x(i)).^2 + (y(i+1)-y(i)).^2).^0.5) ;
            
            % lincoln endurance 2015 and 2016 data has discontinuities at
            % start/finish of track. Ensure smooth transition between start
            % and finish
            if strcmp( obj.event, 'endurance')
                [ x, y ] = obj.smoothEnduranceTrack( ds, x, y, 10 ) ;
                i = 1: length(x)-1 ;
                ds = (((x(i+1)-x(i)).^2 + (y(i+1)-y(i)).^2).^0.5) ;
            end
            
            % get distance of each point from start
            d = zeros( length(x),1) ;
            for i = 2 : length(x)
                d(i) = d(i-1) + ds(i-1) ;
            end
            
            % evenly space points using given track step size
            x = spline( d, x, ( 0 : obj.trackStepSize : max(d) ) ) ;
            y = spline( d, y, ( 0 : obj.trackStepSize : max(d) ) ) ;
            ds = obj.trackStepSize ;
            d = 0 : obj.trackStepSize : max(d) ;
            
            dx = [] ;
            dy = [] ;
            d2x = [] ;
            d2y = [] ;
            % get centered 1st and 2nd derivative dx/d and dy/d
            i = 2 : length(x) - 1 ;
            dx(i) = (x(i+1) - x(i-1)) ./ ( 2*ds) ;
            d2x(i) = (x(i+1) - 2*x(i) + x(i-1)) ./ ds^2 ;
            d2x(end+1) = 0 ;
            dy(i) = (y(i+1) - y(i-1)) ./ ( 2*ds ) ;
            d2y(i) = (y(i+1) - 2*y(i) + y(i-1)) ./ ds^2 ;
            d2y(end+1) = 0 ;
            % take forward derivative at start and backwards derivative at end
            % assume 2nd derivatives = 0 at start and end
            dx(end+1) = (x(end) - x(end-1)) ./ ( ds ) ;
            dx(1) = (x(2) - x(1)) ./ ( ds ) ;
            dy(end+1) = (y(end) - y(end-1)) ./ ( ds ) ;
            dy(1) = (y(2) - y(1)) ./ ( ds ) ;
            
            % get curvature
            k = ( ( dx.*d2y - dy.*d2x ) ./ ( dx.^2 + dy.^2 ).^(3/2) ) ;
            figure
            hold on
            plot(k)
            % filter curvature to smooth data
            [b,a] = butter( 2, obj.butterWorthFilterParameter ) ;
            k = filtfilt(b,a,k) ;
            plot(k)
            
            % create track structure
            track.x = spline( d, x, 0 : obj.trackStepSizeOut : d(end) ) ;
            track.y = spline( d, y, 0 : obj.trackStepSizeOut : d(end) ) ;
            track.d = 0 : obj.trackStepSizeOut : d(end) ;
            track.ds = obj.trackStepSizeOut ;
            track.k = spline( d, k, 0 : obj.trackStepSizeOut : d(end) ) ;
            track.radius = 1./ track.k ;
            track.event = obj.event ;
            if strcmp( obj.event, 'endurance') 
                if ~isempty(obj.numberLaps)
                    track.numberLaps = obj.numberLaps ;
                else
                    track.numberLaps = round( 22000 / max(track.d) ) ;
                end
            end
            if ~strcmp( obj.units, '' )
                track.units = obj.units ;
            end
        end
        
        function [ x, y ] = smoothEnduranceTrack( ~, ds, x, y, tol ) 
            % for endurance track ensure smooth transition between start and end of
            % track
            % at end/start of endurance tracks there are discontinuities,
            % which can be identified due to much larger values of dx
            
            % get median value of track step size
            dsMed = median(ds) ;
            % disregard points with deviation from median greater than tol
            ind = 200 * abs(ds - dsMed) ./ ( ds + dsMed ) < tol ;
            % take end point and start point of remaining points and
            % interpolate linearly between them to reconstruct track
            x = x(ind) ;
            y = y(ind) ;
            
            dist = ( (x(1) - x(end))^2 + (y(1) - y(end))^2 )^0.5 ;
            d = dsMed/dist : dsMed/dist : 1 ;
            
            xi = (x(1)-x(end))*d + x(end) ;
            yi = (y(1)-y(end))*d + y(end) ;
            
            x = [ x, xi ] ;
            y = [ y, yi ] ;
            
        end
    end
    
end

