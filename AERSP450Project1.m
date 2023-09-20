clc ;
clear ;
close all ;

% ---- INPUT VALUES BELOW ---- %

e = 0.5 ;       % eccentricity
T = 8 ;         % orbital period in hours

% ---- INPUT VALUES ABOVE ---- %

% target location (state college)
longitude = -78 ;    % in degrees east negative for west
latitude = 41 ;      % in degrees north negative for south

% constants
lonAN = 60 ;         % cap-omega in degrees 
mu = 3.986*10^5 ;    % mu of earth

% calculations
thetaPeriapsis = 0 ;     % angle theta at Periapsis

T = T*3600 ;    % convert T to seconds

a = (mu*(T/(2*pi))^2)^(1/3) ;   % finding semi-major axis from desired period
p = a*(1-e^2) ;                  % finding  semi-latus rectum from semi-major axis

rp = p/(1+e*cos(thetaPeriapsis)) ;                  % radius of periapsis for designated conditions

C1 = tand(longitude) ;
C2 = tand(latitude)*sqrt(1+C1^2) ;

R_ECI(1) = rp/sqrt(1+C1^2+C2^2) ;
R_ECI(2) = R_ECI(1)*tand(longitude)  ;
R_ECI(3) = tand(latitude)*sqrt((R_ECI(1)^2+R_ECI(2)^2)) 

V_Perifocal(1) = sqrt(mu/p)*-sin(thetaPeriapsis) ;
V_Perifocal(2) = sqrt(mu/p)*(e+cos(thetaPeriapsis)) ;
V_Perifocal(3) = 0 

R_Perifocal(1) = rp*cos(thetaPeriapsis) ;
R_Perifocal(2) = rp*sin(thetaPeriapsis) ;
R_Perifocal(3) = 0 

h = cross(R_Perifocal,V_Perifocal)
