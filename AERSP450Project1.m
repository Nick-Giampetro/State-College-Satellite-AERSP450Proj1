

mu = 3.986*10^5 ;    % mu of earth

longitude = -78 ;    % in degrees east negative for west
latitude = 41 ;      % in degrees north negative for south

T = 8 ;         % orbital period in hours
T = T*3600 ;    % convert T to seconds

e = 0.5 ;       % e

a = (mu*(T/(2*pi))^2)^(1/3)

rp = a*(1-e)
