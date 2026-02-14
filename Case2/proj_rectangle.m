function y_proj = proj_rectangle(theta, y, theta_min, theta_max, delta)
% 向量化实现的光滑投影算子
% 输入:
%   theta: 当前参数向量 (n x 1)
%   y: 原始更新方向向量 (n x 1)
%   theta_min: 参数下界向量 (n x 1)
%   theta_max: 参数上界向量 (n x 1)
%   delta: 边界过渡区厚度 (标量)
% 输出:
%   y_proj: 投影后的更新方向向量 (n x 1)

% 验证输入维度一致性
% if ~isequal(size(theta), size(y), size(theta_min), size(theta_max))
%     error('所有输入向量必须具有相同维度');
% end

% 检查delta有效性
% if delta <= 0
%     error('delta必须为正数');
% end

% 验证参数范围有效性
% range_check = (theta_max - theta_min) < 2 * delta;
% if any(range_check)
%     error('每个参数的范围必须大于2*delta');
% end

% 创建逻辑索引
lower_bound = (theta <= theta_min + delta) & (y < 0);  % 下边界过渡区
upper_bound = (theta >= theta_max - delta) & (y > 0);  % 上边界过渡区
safe_region = ~(lower_bound | upper_bound);            % 安全区域

% 初始化输出向量
y_proj = zeros(size(y));

% 下边界投影计算
y_proj(lower_bound) = y(lower_bound) .* ...
    (theta(lower_bound) - theta_min) / delta;

% 上边界投影计算
y_proj(upper_bound) = y(upper_bound) .* ...
    (theta_max - theta(upper_bound)) / delta;

% 安全区域保持原值
y_proj(safe_region) = y(safe_region);
end