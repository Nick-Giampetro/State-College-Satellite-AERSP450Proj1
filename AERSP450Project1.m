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

% calculation to get R peripasis in ECEF
C1 = tand(longitude) ;
C2 = tand(latitude)*sqrt(1+C1^2) ;
R_ECEF(1) = rp/sqrt(1+C1^2+C2^2) ;
R_ECEF(2) = R_ECEF(1)*tand(longitude)  ;
R_ECEF(3) = tand(latitude)*sqrt((R_ECEF(1)^2+R_ECEF(2)^2)) 

% Perifocal coordinate calculations
V_Perifocal(1) = sqrt(mu/p)*-sind(thetaPeriapsis) ;
V_Perifocal(2) = sqrt(mu/p)*(e+cosd(thetaPeriapsis)) ;
V_Perifocal(3) = 0 ;
 
R_Perifocal(1) = rp*cosd(thetaPeriapsis) ;
R_Perifocal(2) = rp*sind(thetaPeriapsis) ;
R_Perifocal(3) = 0 ;

% argument of assending node calculation
w = asind(R_ECEF(3)/(R_Perifocal(1)*sind(I)))

% GMST calculation
C1 = cosd(lonAN)/R_ECEF(2) + sind(lonAN)/R_ECEF(1) ;
C2 = tand(lonAN)*cosd(lonAN)/R_ECEF(2) - sind(lonAN)/(R_ECEF(1)*tand(lonAN)) ;
C3 = R_ECEF(1)/R_ECEF(2) + R_ECEF(2)/R_ECEF(1) ;
GMST = acosd((R_Perifocal(1) * (cosd(w)*C1 - cosd(I)*sind(w)*C2))/C3) 

% ECI vector transformation
cEP = dcm3axis(w)*dcm1axis(I)*dcm3axis(lonAN) ;
cPE = cEP';
R_ECI = cPE * R_Perifocal' ;
V_ECI = cPE * V_Perifocal' ;

% ODE45 and FG orbit propagation
init = [ R_ECI(1) R_ECI(2) R_ECI(3) V_ECI(1) V_ECI(2) V_ECI(3)] ;
totalT = T*24 ;
t = linspace(1,totalT,25000) ;
options = odeset('reltol',1e-12,'abstol',1e-12);
[t,R_ODE45] = ode45( @(t,R_ODE45) TwoBP(t,R_ODE45,mu) , t , init, options) ;
[R_FG, V_FG] = FGProp(t,R_ECI,V_ECI,e,a,p,mu) ;

% error calculations
r1 = zeros(length(t),3) ;
v1 = zeros(length(t),3) ;
for idx = 1:3
        r1(:,idx) = R_ODE45(:,idx) ;
        v1(:,idx) = R_ODE45(:,idx+3) ;
end
R_Error = abs(R_FG - r1);
V_Error = abs(V_FG - v1) ;

% ECEF coordinate transformation
omega = 360/24/3600 ;
R_ODE45_ECEF = ECI2ECEF(r1, omega, t, GMST) ;
R_FG_ECEF = ECI2ECEF(R_FG, omega, t, GMST) ;

% ground track plotting
grdTrck_ODE45 = groundTrack(t,R_ODE45,GMST) ;
grdTrck_FG = groundTrack(t,R_FG,GMST) ;

% finding orbital elements at each time
orbital_elements = zeros(length(t),6) ;
for idx = 1:length(t)    
    [orbital_elements(idx,1), orbital_elements(idx,2), orbital_elements(idx,3), orbital_elements(idx,4), orbital_elements(idx,5), orbital_elements(idx,6)] = RV2OE(r1(idx,:),v1(idx,:),mu) ;
end

% graphs
f = figure ;
subplot(1,1,1)
plot3(R_ODE45_ECEF(:,1), R_ODE45_ECEF(:,2), R_ODE45_ECEF(:,3), 'r', R_FG_ECEF(:,1), R_FG_ECEF(:,2), R_FG_ECEF(:,3), 'b')
xlabel('X (KM)')
ylabel('Y (KM)')
zlabel('Z (KM)')
exportgraphics(f,['3D' '.jpg'])

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



% propogates F and G functions to plot orbit
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


% calculates r and v vectors from F and G functions
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


% calculates ground track plot
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


% Two body problem function
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


% ECI to ECEF translation
function r_ecef = ECI2ECEF(r_eci, omega, t, GMST)
    r_ecef = zeros(length(t),3) ;

    for idx = 1 : length(t)
        gamma = GMST + t(idx) * omega ;
        cECEF = dcm3axis(gamma) ;
        r_ecef(idx,:) = cECEF' * r_eci(idx,:)' ;
    end
end


% Newton Raphson method
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


% gets orbital elemets from an R and V value
function [a,e,I,RAAN,AOP,f] = RV2OE(r,v,mu)
        
    % declaring gravitational constant of earth and time provided and ECI
    ECI = [[1 0 0];[0 1 0];[0 0 1]];
    
    R1 = norm(r);
    V1 = norm(v);

    energy = (V1^2)/2-mu/R1;
    
    h = cross(r,v);
    H = norm(h);
    
    p = H^2/mu;
    
    % calculating semi major axis
    a = -mu/(2*energy);
    
    % calculating eccentricity 
    eV = cross(v,h)/mu - (r/R1);
    e = norm(eV);
    
    % calculating orbital inclination
    I = acos(dot(ECI(3,:),h)/H);
    
    % calculating longitude of the ascending node
    n = cross(ECI(3,:),h)/norm(cross(ECI(3,:),h));
    
    if(dot(n,ECI(1,:)) >= 0)
        RAAN = atan(dot(n,ECI(2,:))/dot(n,ECI(1,:)));
    elseif (dot(n,ECI(1,:)) < 0)
        RAAN = atan(dot(n,ECI(2,:))/dot(n,ECI(1,:)))+pi;
    end
    
    % calculating argument of periapsis
    if(dot(eV,ECI(:,3)) >= 0)
        AOP = acos(dot(eV,n)/e);
    elseif (dot(eV,ECI(:,3)) < 0)
        AOP = -acos(dot(eV,n)/e);
    end

    if(dot(r,eV) >= 0)
        f = acosd(dot(eV,r)/(e*R1)) ;
    elseif(dot(r,eV) < 0)
        f = 360 - acosd(dot(eV,r)/(e*R1)) ;
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
