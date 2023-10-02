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

% calculation of constants needed to find ECEF position vector
C1 = tand(longitude) ;
C2 = tand(latitude)*sqrt(1+C1^2) ;

% calculation of position vector for question 3 of analytical part using C1
% and C2 from above along with latitude and longitude of State College
R_ECEF(1) = rp/sqrt(1+C1^2+C2^2) ;
R_ECEF(2) = R_ECEF(1)*tand(longitude)  ;
R_ECEF(3) = tand(latitude)*sqrt((R_ECEF(1)^2+R_ECEF(2)^2))  

V_Perifocal(1) = sqrt(mu/p)*-sind(thetaPeriapsis) ;
V_Perifocal(2) = sqrt(mu/p)*(e+cosd(thetaPeriapsis)) ;
V_Perifocal(3) = 0 ;  
 
% calculation of position vector in orbital frame needed for question 4b of
% analytical part
R_Perifocal(1) = rp*cosd(thetaPeriapsis) ;
R_Perifocal(2) = rp*sind(thetaPeriapsis) ;
R_Perifocal(3) = 0 ;    

% calculation of argument of periapsis needed to find GMST using sinw
% expression found in part 4c of analytical part of project
w = asind(R_ECEF(3)/(R_Perifocal(1)*sind(I)))

% constants that can be used in GMST calculation in order to prevent
% complexity and potential syntax mistakes in GMST calculation
C1 = cosd(lonAN)/R_ECEF(2) + sind(lonAN)/R_ECEF(1) ;
C2 = tand(lonAN)*cosd(lonAN)/R_ECEF(2) - sind(lonAN)/(R_ECEF(1)*tand(lonAN)) ;
C3 = R_ECEF(1)/R_ECEF(2) + R_ECEF(2)/R_ECEF(1) ;

% GMST calculation needed for question 4d of analytical part of project
GMST = acosd((R_Perifocal(1) * (cosd(w)*C1 - cosd(I)*sind(w)*C2))/C3) 

% calculation of rotational matrix needed to find inertial vectors and the
% calculations of the inertial position and velocity vectors needed for
% questions 4f and 5 of analytical part
cEP = dcm3axis(w)*dcm1axis(I)*dcm3axis(lonAN) ;
cPE = cEP';
R_ECI = cPE * R_Perifocal' ;
V_ECI = cPE * V_Perifocal' ;

% application of ODE45 for propogating orbit for question one of numerical
% part of project
init = [ R_ECI(1) R_ECI(2) R_ECI(3) V_ECI(1) V_ECI(2) V_ECI(3)] ;
totalT = T*24 ;
t = linspace(1,totalT,25000) ;
options = odeset('reltol',1e-12,'abstol',1e-12);
[t,R_ODE45] = ode45( @(t,R_ODE45) TwoBP(t,R_ODE45,mu) , t , init, options) ;

% application of FG propgation function needed for question 4 of numerical
% part of project
[R_FG, V_FG] = FGProp(t,R_ECI,V_ECI,e,a,p,mu) ;


r1 = zeros(length(t),3) ;
v1 = zeros(length(t),3) ;
for idx = 1:3
        r1(:,idx) = R_ODE45(:,idx) ;
        v1(:,idx) = R_ODE45(:,idx+3) ;
end

% anular velocity of the earth calculation
omega = 360/24/3600 ;

R_ODE45_ECEF = ECI2ECEF(r1, omega, t, GMST) ;
R_FG_ECEF = ECI2ECEF(R_FG, omega, t, GMST) ;

% calculation of errors between two methods needed for question 5 of
% numerical part of project
R_Error = abs(R_FG - r1);
V_Error = abs(V_FG - v1) ;

% application of groundtrack function created to answer question 3 of the
% numerical part of project
grdTrck_ODE45 = groundTrack(t,R_ODE45,GMST) ;
grdTrck_FG = groundTrack(t,R_FG,GMST) ;

% 3D plot of ecef framed orbit needed for question 2a of numerical part
f = figure ;
subplot(1,1,1)
plot3(R_ODE45_ECEF(:,1), R_ODE45_ECEF(:,2), R_ODE45_ECEF(:,3), 'r', R_FG_ECEF(:,1), R_FG_ECEF(:,2), R_FG_ECEF(:,3), 'b')
xlabel('X (KM)')
ylabel('Y (KM)')
zlabel('Z (KM)')
exportgraphics(f,['3D' '.jpg'])

% figures for questions 2b, 2c, and 2d of numerical part
f = figure ;
subplot(1,1,1)
plot(R_ODE45_ECEF(:,1), R_ODE45_ECEF(:,2), 'r', R_FG_ECEF(:,1), R_FG_ECEF(:,2), 'b')
xlabel('X (KM)')
ylabel('Y (KM)')
exportgraphics(f,['X vs Y' '.jpg'])
f = figure ;
subplot(1,1,1)
plot(R_ODE45_ECEF(:,2), R_ODE45_ECEF(:,3), 'r', R_FG_ECEF(:,2), R_FG_ECEF(:,3), 'b')
xlabel('Y (KM)')
ylabel('Z (KM)')
exportgraphics(f,['Y vs Z' '.jpg'])

f = figure ;
subplot(1,1,1)
plot(R_ODE45_ECEF(:,1), R_ODE45_ECEF(:,3), 'r', R_FG_ECEF(:,1), R_FG_ECEF(:,3), 'b')
xlabel('X (KM)')
ylabel('Z (KM)')
exportgraphics(f,['X vs Z' '.jpg'])

% imported coastlines of earth for groundtrack plot and plotted
% groundtracks needed for question 3 of numerical part of project
f = figure ;
load('topo.mat','topo');
topoplot = [topo(:,181:360),topo(:,1:180)];
contour(-180:179,-90:89,topoplot,[0,0],'black','linewidth',1);
grid on
grid minor
axis equal
hold on
plot( grdTrck_ODE45(:,1), grdTrck_ODE45(:,2),'b', grdTrck_FG(:,1), grdTrck_FG(:,2), 'g', longitude, latitude,'X')% <- ENTER THE LONGITUDE AND LATITUDE YOU HAVE COMPUTED HERE
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
set(gca,'FontSize',18)
exportgraphics(f,['long vs lat' '.jpg'])



% function for F and G function propagation 
function [r,v] = FGProp(t,r0,v0,e,a,p,mu)
    
    r = zeros(length(t),3) ;
    v = zeros(length(t),3) ;
    
    r(1,:) = r0 ;
    r(2,:) = v0 ;

    E0 = NewtonMethod(a,e,mu,t(1),0,0.000001) ;
    
    for idx = 2:length(t) 
        E = NewtonMethod(a,e,mu,t(idx),0,0.0000001) ;
    
        dE = E - E0 ;

        [r(idx,:),v(idx,:)] = FGFunc(r0,v0,dE,t(idx),e,p,a,mu) ;

    end
end



function  [r,v] = FGFunc(r0,v0,dE,dt,e,p,a,mu)
R0 = norm(r0) ;

theta0 = acos((1/e)*(p/R0-1)) ;
R = p/(1+e*cos(theta0+dE)) ;

f = 1 - a/R0 * (1 - cos(dE)) ;
g = dt - sqrt(a^3/mu) * (dE-sin(dE)) ;
fdot = (-1*sqrt(mu*a)*sin(dE)) / (R*R0) ;
gdot = 1 - a/R *(1-cos(dE)) ;

r = f*r0 + g*v0 ;
v = fdot*r0 + gdot*v0 ;
end


% function to calulate ground tracks of satellite
function gt = groundTrack(t,r,GMST)
    r1 = zeros(length(t),3) ;

    for idx = 1:3
        r1(:,idx) = r(:,idx) ;
    end

    e_omega = 360/24/3600 ;

    r_ecef = ECI2ECEF(r1,e_omega,t,GMST) ;
    
    gt = zeros(length(t),2) ;
    phase = 0 ;

    for idx = 1:length(t)
        if mod(phase,3) == 0
            gt(idx,1) = atand(r_ecef(idx,2)/r_ecef(idx,1)) ;
            gt(idx,2) = atand( r_ecef(idx,3) / sqrt( r_ecef(idx,1)^2 + r_ecef(idx,2)^2 )) ;
        elseif  mod(phase,3) == 1
            gt(idx,1) = 180 + atand(r_ecef(idx,2)/r_ecef(idx,1)) ;
            gt(idx,2) = atand(r_ecef(idx,3)/sqrt(r_ecef(idx,1)^2+r_ecef(idx,2)^2)) ;
        elseif  mod(phase,3) == 2
            gt(idx,1) = -180 + atand(r_ecef(idx,2)/r_ecef(idx,1)) ;
            gt(idx,2) = atand(r_ecef(idx,3)/sqrt(r_ecef(idx,1)^2+r_ecef(idx,2)^2)) ;
        end
           
        if (idx > 1) && (gt(idx,1)-gt(idx-1,1) < -90 ) && (gt(idx-1) > 0)
           phase = phase + 1 ; 
           gt(idx,1) = 180 + atand(r_ecef(idx,2)/r_ecef(idx,1)) ;
        elseif (idx > 1) && (gt(idx,1) > 180)
           phase = phase + 1 ;
           gt(idx,1) = -180 + atand(r_ecef(idx,2)/r_ecef(idx,1)) ;
        elseif (idx > 1) && (gt(idx,1)-gt(idx-1,1) < -90) && (gt(idx-1) < 0)
           phase = phase + 1 ; 
           gt(idx,1) = atand(r_ecef(idx,2)/r_ecef(idx,1)) ;
        end
    end
end


% function to code eqns of motion for two body program so ode45 can be used
function dx = TwoBP(~,r,mu)
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


% function to calculate eci vectors from ecef framed vectors
function r_ecef = ECI2ECEF(r_eci, omega, t, GMST)
    r_ecef = zeros(length(t),3) ;

    for idx = 1 : length(t)
        gamma = GMST + t(idx) * omega ;
        cECEF = dcm3axis(gamma) ;
        r_ecef(idx,:) = cECEF' * r_eci(idx,:)' ;
    end
end


% function to perform the NewtonRaphson iteration for F&G propagation
function an = NewtonMethod(a,e,mu,t,tp,delta)

M = sqrt(mu/(a^3))*(t-tp) ;
error = delta + 1 ;
Eo = M ;

while error > delta
    En = Eo - ((M+e*sin(Eo)-Eo)/(e*cos(Eo)-1)) ;
    error = abs(En-Eo) ;
    Eo = En ;
end

an = En;

end


% creates a dcm for an angle about axis 1
function r = dcm1axis(ang)
r = [1 0 0 ; 0 cosd(ang) sind(ang) ; 0 -sind(ang) cosd(ang)];
end

% creates a dcm for an angle about axis 3
function r = dcm3axis(ang)
r = [cosd(ang) sind(ang) 0 ; -sind(ang) cosd(ang) 0 ; 0 0 1];
end
