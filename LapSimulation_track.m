 
%% Load Track .mat File

load('Lincoln2016_Endurance.mat')
% load('track_FSG_2013.mat')
% for x = 1:(length(GPS_Latitude.Value) - 3)
x = 1:(length(GPS_Latitude.Value) - 3);

    % Get lengths between points
a = sqrt( (GPS_Latitude.Value(x+1) - GPS_Latitude.Value(x)).^2 + (GPS_Longitude.Value(x+1) - GPS_Longitude.Value(x)).^2); 
b = sqrt( (GPS_Latitude.Value(x+2) - GPS_Latitude.Value(x+1)).^2 + (GPS_Longitude.Value(x+2) - GPS_Longitude.Value(x+1)).^2);
c = sqrt( (GPS_Latitude.Value(x+3) - GPS_Latitude.Value(x+2)).^2 + (GPS_Longitude.Value(x+3) - GPS_Longitude.Value(x+2)).^2);

    % Get Angles
A = acos( (b(x).^2 + c(x).^2 + a(x).^2)./(2*b(x).*c(x)));

    % Get radii
Radius = a(x)./(2*sin(180 - A(x))); 

% end

%% Build Working Track Variables 

   % track.dx        = 0.75;
   % track.xbase     = 0:track.dx:max(track.x);
   % track.crv       = interp1(track.x,track.crv_filt,track.xbase,'spline');

   % track.x         = track.xbase;
   
  % SectorLength = zeros(1:length(GPS_Latitude.Value)-1);
   
    
%  for x = 1:length(GPS_Latitude.Value - 1)
 %     
  %    SectorLength(x) = sqrt((GPS_Latitude.Value(x+1) - GPS_Latitude.Value(x))^2 + (GPS_Longitude.Value(x+1) - GPS_Longitude.Value(x))^2);
   %   Angle(x) = (tan((GPS_Longitude.Value(x+1) - GPS_Longitude.Value(x))/(GPS_Latitude.Value(x+1) - GPS_Latitude.Value(x))));
      
  %end
  
  
      
      
    %% Calculate Instantaneous Radius
    
  %  track.radius    = (1 ./ track.crv);
  %  plot(track.radius);
   % ylim([-500, 500]);