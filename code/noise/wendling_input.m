function [noise, raw] = wendling_input(cfg, node_count, steps)
% WENDLING_INPUT  用另一个 Wendling 模型的输出 y_out 充当外部驱动 p(t)。
%
%   [noise, raw] = wendling_input(cfg, node_count, steps)
%
% 输入
%   cfg         配置结构体（default_config），本函数用到的字段：
%                 noise_wendler_params  驱动模型的 [G_e G_s G_f]，默认 [4.25 22 10]
%                 noise_wendler_dynamic 驱动模型是否启用 ESF 增益动力学，默认 0
%                 noise_wendler_normalize 是否把驱动输出线性映射到噪声量纲，默认 true
%                 noise_mean / noise_sigma / simulation_time / dt  同 gauss_noise
%   node_count  节点数（每个节点一个独立驱动模型实例）
%   steps       需要的样本数
%
% 输出
%   noise       node_count × steps，单位 pulses/s（可直接喂给 simulate_network 的模型）
%   raw         node_count × steps，驱动模型的原始 y_out（mV），未做量纲映射
%
% 做什么
%   这是一个"替代驱动源"的实验分支：不用白噪声 p(t)，而是先用一个 Wendling
%   模型（默认 G_e = 4.25，即已经越过 ESF 阈值的兴奋性增益）自由跑一段，
%   把它产生的 EEG 样信号当作下游节点的传入驱动。
%
%   为什么要做量纲映射
%     y_out 是 mV 量级的电位差（几到几十），而 p(t) 在原文里是脉冲密度
%     （均值 90 pulses/s）。两者直接相接会把下游打到饱和或完全驱动不起来。
%     默认 normalize = true：把 y_out 去均值、按噪声标准差缩放、再加噪声均值，
%     即 p_i(t) = mean + sigma * (y(t) - mean(y)) / std(y)。
%     这样"驱动的时间结构来自模型、驱动强度仍由 noise_mean/sigma 控制"，
%     换驱动源时不必重调下游参数。想保留原始 mV 尺度就把 normalize 设 false。
%
%   用途定位
%     纯粹是实验分支，论文的三张图都不用这条路径（默认 noise_type = 'gauss'）。
%     如果以后要论证"耦合进来的信号本身是振荡性的、与白噪声驱动不同"，
%     这个函数就是现成的对照条件。
%
% 调用关系
%   上游：simulate_network（noise_type = 'wendling'）
%   下游：wendling（驱动模型）、gauss_noise（驱动模型自己的 p(t)）
%
% 备注
%   驱动模型的初始状态由 wendling 构造函数里的 rand 决定，因此本函数同样受
%   simulate_network 开头的 rng(cfg.seed) 控制，可复现。

if nargin < 2 || isempty(node_count), node_count = 1; end

dt = cfg.dt;

params  = getfield_default(cfg, 'noise_wendler_params', [4.25 22 10]);
if_dy   = getfield_default(cfg, 'noise_wendler_dynamic', 0);
norm_on = getfield_default(cfg, 'noise_wendler_normalize', true);
mu      = getfield_default(cfg, 'noise_mean', 90);
sigma   = getfield_default(cfg, 'noise_sigma', 30);

if numel(params) ~= 3
    error('cfg.noise_wendler_params 需要 [G_e G_s G_f] 三个数，当前 %d 个', numel(params));
end

% 驱动模型自己的输入：仍用高斯噪声（分段设置沿用 cfg，一般就是单段同分布）
drive = gauss_noise(cfg, node_count, steps);

raw   = zeros(node_count, steps);

for n = 1:node_count
    m = wendling(params(1), params(2), params(3), dt, if_dy);
    for t = 1:steps
        [~, ~, y_out] = m.simulate(drive(n, t));
        raw(n, t) = y_out;
    end
end

if norm_on
    noise = zeros(node_count, steps);
    for n = 1:node_count
        x  = raw(n, :);
        sx = std(x);
        if sx < eps
            % 驱动模型没跑起来（全常数）：直接退回高斯噪声，避免除零后全 NaN
            noise(n, :) = mu + sigma * randn(1, steps);
        else
            noise(n, :) = mu + sigma * (x - mean(x)) / sx;
        end
    end
else
    noise = raw;
end
end

% ======================================================================
function v = getfield_default(cfg, name, default)
if isfield(cfg, name) && ~isempty(cfg.(name))
    v = cfg.(name);
else
    v = default;
end
end
