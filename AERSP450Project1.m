clc ;
clear ;
close all ;

% ---- INPUT VALUES BELOW ---- %

e = 0.2;        % eccentricity
T = 2 ;         % orbital period in hours
I = 55;         % must be greater than the latitude value below (line 15)

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

a = (mu*(T/(2*pi))^2)^(1/3) ;    % finding semi-major axis from desired period
p = a*(1-e^2) ;                  % finding  semi-latus rectum from semi-major axis

rp = p/(1+e*cos(thetaPeriapsis)) ;                  % radius of periapsis for designated conditions

C1 = tand(longitude) ;
C2 = tand(latitude)*sqrt(1+C1^2) ;

R_ECEF(1) = rp/sqrt(1+C1^2+C2^2) ;
R_ECEF(2) = R_ECEF(1)*tand(longitude)  ;
R_ECEF(3) = tand(latitude)*sqrt((R_ECEF(1)^2+R_ECEF(2)^2)) 

% V_Perifocal(1) = sqrt(mu/p)*-sind(thetaPeriapsis) ;
% V_Perifocal(2) = sqrt(mu/p)*(e+cosd(thetaPeriapsis)) ;
% V_Perifocal(3) = 0 ;
% 
R_Perifocal(1) = rp*cosd(thetaPeriapsis) ;
R_Perifocal(2) = rp*sind(thetaPeriapsis) ;
R_Perifocal(3) = 0 

w = asind(R_ECEF(3)/(R_Perifocal(1)*sind(I)))

C1 = cosd(lonAN)/R_ECEF(2) + sind(lonAN)/R_ECEF(1) ;
C2 = tand(lonAN)*cosd(lonAN)/R_ECEF(2) - sind(lonAN)/(R_ECEF(1)*tand(lonAN)) ;
C3 = R_ECEF(1)/R_ECEF(2) + R_ECEF(2)/R_ECEF(1) ;

GMST = acosd((R_Perifocal(1) * (cosd(w)*C1 - cosd(I)*sind(w)*C2))/C3)


% cEP = dcm3axis(w)*dcm1axis(I)*dcm3axis(lonAN) ;
% cPE = cEP';
% R_ECI = cPE * R_Perifocal'








% creates a dcm for an angle about axis 1
function r = dcm1axis(ang)
r = [1 0 0 ; 0 cosd(ang) sind(ang) ; 0 -sind(ang) cosd(ang)];
end

% creates a dcm for an angle about axis 3
function r = dcm3axis(ang)
r = [cosd(ang) sind(ang) 0 ; -sind(ang) cosd(ang) 0 ; 0 0 1];
end

