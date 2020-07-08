function g_n = gravity(lat, h)
% gravity: calculates gravity vector in the navigation frame.
%
% INPUT:
%       lat: Mx1 latitude (radians).
%         h: Mx1 altitude (m).
%
% OUTPUT:
%		g_n: Mx1 gravity vector in the nav-frame (m/s^2).

h = abs(h);
sin1 = sin(lat);
sin2 = sin(2.*lat);

g0 = 9.780318 * ( 1 + 5.3024e-03.*(sin1).^2 - 5.9e-06.*(sin2).^2 );

[RM,RN] = radius(lat);

Ro = sqrt(RN .* RM);

g = (g0 ./ (1 + (h ./ Ro)).^2);

Z = zeros(size(lat));

g_n = [Z Z g];

end

