function data_r = p2r(data_p, r, phi, n)
% function data_r = p2r(data_p, r, phi, n)
%
% Convert polar coordinate data to cartisian coordinate data.
% phi is in degrees.

INNER_VALUE = min(data_p(:));
OUTER_VALUE = INNER_VALUE;

data_r = zeros(n);

center = (n-1)/2;

for i=0:(n-1)
  for j=0:(n-1)
    x = (j - center)/n * r(end) + r(end)/(2*n);
    y = (i - center)/n * r(end) + r(end)/(2*n);
    
    radius = sqrt(x^2 + y^2);
    angle = atan2(-y, x) * 180 / pi;
    
    if (radius < r(1))
      data_r(i+1, j+1) = INNER_VALUE;
    elseif (radius > r(end))
      data_r(i+1, j+1) = OUTER_VALUE;
    else
      I = find(r <= radius);
      J = find(phi > angle);
      
      if (length(I) < 1 | length(J) < 1)
        radius
        angle
        error('Assertion failed');
      end
      
      data_r(i+1, j+1) = data_p(I(end), J(1));
    end
  end
end