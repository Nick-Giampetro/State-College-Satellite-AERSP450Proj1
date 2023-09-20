

mu = 3.986*10^5 ;    % mu of earth

longitude = -78 ;    % in degrees east negative for west
latitude = 41 ;      % in degrees north negative for south

T = 8 ;         % orbital period in hours
T = T*3600 ;    % convert T to seconds

e = 0.5 ;       % e

a = (mu*(T/(2*pi))^2)^(1/3)

rp = a*(1-e)

C1 = tand(longitude) ;
C2 = tand(latitude)*sqrt(1+C1^2) ;

R(1) = rp/sqrt(1+C1^2+C2^2);
R(2) = R(1)*tand(longitude) ;
R(3) = tand(latitude)*sqrt((R(1)^2+R(2)^2))
