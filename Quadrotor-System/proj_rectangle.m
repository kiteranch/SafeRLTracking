function y_proj = proj_rectangle(theta, y, theta_min, theta_max, delta)
% Vectorized smooth projection operator
% Input:
%   theta: Current parameter vector (n x 1)
%   y: Original update direction vector (n x 1)
%   theta_min: Parameter lower bound vector (n x 1)
%   theta_max: Parameter upper bound vector (n x 1)
%   delta: Boundary transition zone thickness (scalar)
% Output:
%   y_proj: Projected update direction vector (n x 1)

% Verify input dimension consistency
% if ~isequal(size(theta), size(y), size(theta_min), size(theta_max))
%     error('All input vectors must have the same dimension');
% end

% Check delta validity
% if delta <= 0
%     error('delta must be positive');
% end

% Verify parameter range validity
% range_check = (theta_max - theta_min) < 2 * delta;
% if any(range_check)
%     error('The range of each parameter must be greater than 2*delta');
% end

% Create logical indices
lower_bound = (theta <= theta_min + delta) & (y < 0);  % Lower boundary transition zone
upper_bound = (theta >= theta_max - delta) & (y > 0);  % Upper boundary transition zone
safe_region = ~(lower_bound | upper_bound);            % Safe region

% Initialize output vector
y_proj = zeros(size(y));

% Lower boundary projection calculation
y_proj(lower_bound) = y(lower_bound) .* ...
    (theta(lower_bound) - theta_min) / delta;

% Upper boundary projection calculation
y_proj(upper_bound) = y(upper_bound) .* ...
    (theta_max - theta(upper_bound)) / delta;

% Safe region keep original value
y_proj(safe_region) = y(safe_region);
end