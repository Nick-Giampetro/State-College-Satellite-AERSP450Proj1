clc ;
clear ;
close all ;

% ---- INPUT VALUES BELOW ---- %

e = 0.1;        % eccentricity
T = 2 ;         % orbital period in hours
I = 45;         % must be greater than the latitude value below (line 15)

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

V_Perifocal(1) = sqrt(mu/p)*-sind(thetaPeriapsis) ;
V_Perifocal(2) = sqrt(mu/p)*(e+cosd(thetaPeriapsis)) ;
V_Perifocal(3) = 0 ;
 
R_Perifocal(1) = rp*cosd(thetaPeriapsis) ;
R_Perifocal(2) = rp*sind(thetaPeriapsis) ;
R_Perifocal(3) = 0 

w = asind(R_ECEF(3)/(R_Perifocal(1)*sind(I)))

C1 = cosd(lonAN)/R_ECEF(2) + sind(lonAN)/R_ECEF(1) ;
C2 = tand(lonAN)*cosd(lonAN)/R_ECEF(2) - sind(lonAN)/(R_ECEF(1)*tand(lonAN)) ;
C3 = R_ECEF(1)/R_ECEF(2) + R_ECEF(2)/R_ECEF(1) ;

GMST = acosd((R_Perifocal(1) * (cosd(w)*C1 - cosd(I)*sind(w)*C2))/C3)


cEP = dcm3axis(w)*dcm1axis(I)*dcm3axis(lonAN) ;
cPE = cEP';
R_ECI = cPE * R_Perifocal' ;
V_ECI = cPE * V_Perifocal' ;

init = [ R_ECI(1) R_ECI(2) R_ECI(3) V_ECI(1) V_ECI(2) V_ECI(3)] ;
totalT = T*24 ;
t = linspace(1,totalT,100000) ;
options = odeset('reltol',1e-12,'abstol',1e-12);

[t,R_ODE45] = ode45( @(t,R_ODE45) TwoBP(t,R_ODE45,mu) , t , init, options) ;

grdTrck = groundTrack(t,R_ODE45,GMST,mu) ;

f = figure ;
subplot(1,1,1)
plot(R_ODE45(:,1),R_ODE45(:,2),'r')
exportgraphics(f,['X vs Y' '.jpg'])
f = figure ;
subplot(1,1,1)
plot(R_ODE45(:,2),R_ODE45(:,3),'g')
exportgraphics(f,['Y vs Z' '.jpg'])
f = figure ;
subplot(1,1,1)
plot(R_ODE45(:,1),R_ODE45(:,3),'b')
exportgraphics(f,['X vs Z' '.jpg'])

f = figure ;
subplot(1,1,1)
plot(grdTrck(:,1),grdTrck(:,2),'b')
exportgraphics(f,['long vs lat' '.jpg'])


function dx = TwoBP(t,r,mu)
    x = r(1) ;
    y = r(2) ;
    z = r(3) ;
    xdot = r(4) ;
    ydot = r(5) ;
    zdot = r(6) ;

    R = sqrt(x^2 + y^2 + z^2) ;

    xdoubledot = -mu/R^3 * x ;
    ydoubledot = -mu/R^3 * y ;
    zdoubledot = -mu/R^3 * z ;

    dx = [ xdot ; ydot ; zdot ; xdoubledot ; ydoubledot ; zdoubledot ] ;
end


function gt = groundTrack(t,r,GMST,mu)
    
    for idx = 1:3
        r1(:,idx) = r(:,idx) ;
        v1(:,idx) = r(:,idx+3) ;
    end

    %e_omega = 360/24/3600 ;
    e_omega=0;

    [r_ecef,v_ecef] = ECI2ECEF(r1,v1,e_omega,t,GMST,mu) ;
    
    gt = zeros(length(t),2) ;
    crt = 0 ;

    for idx = 1:length(t)
        gt(idx,1) = crt + atand(r_ecef(idx,2)/r_ecef(idx,1)) ;
        gt(idx,2) = atand(r_ecef(idx,3)/sqrt(r_ecef(idx,1)^2+r_ecef(idx,2)^2)) ;
        if (idx > 1) && (gt(idx,1)-gt(idx-1,1) > 90)
            crt = 180 ;   
        end
    end
end

function [r_ecef,v_ecef] = ECI2ECEF(r_eci, v_eci, omega, t, GMST, mu)
    r_ecef = zeros(length(t),3) ;
    v_ecef = zeros(length(t),3) ;

    for idx = 1 : length(t)
        gamma = GMST + t(idx) * omega ;
        cECEF = dcm1axis(gamma) ;
        r_ecef(idx,:) = cECEF * r_eci(idx,:)' ;
        v_ecef(idx,:) = cECEF * v_eci(idx,:)' ;
    end
end




% creates a dcm for an angle about axis 1
function r = dcm1axis(ang)
r = [1 0 0 ; 0 cosd(ang) sind(ang) ; 0 -sind(ang) cosd(ang)];
end

% creates a dcm for an angle about axis 3
function r = dcm3axis(ang)
r = [cosd(ang) sind(ang) 0 ; -sind(ang) cosd(ang) 0 ; 0 0 1];
end

