function noise = gauss_noise(cfg, node_count, steps)
% GAUSS_NOISE  生成 N 个节点的外部驱动 p(t)：分段高斯噪声（同一节点内首尾拼接）。
%
%   noise = gauss_noise(cfg, node_count, steps)
%
% 输入
%   cfg         配置结构体（default_config），本函数用到的字段：
%                 simulation_time   总时长 T（s）
%                 dt                步长（s）
%                 noise_mean        默认均值（pulses/s）
%                 noise_sigma       默认标准差
%                 noise_segments    各段占总时长的比例，默认 [0.3 0.5 0.2]
%                 noise_seg_means   （可选）逐段均值，缺省 = 全部用 noise_mean
%                 noise_seg_sigmas  （可选）逐段标准差，缺省 = 全部用 noise_sigma
%                 noise_method      'normal' | 'truncated' | 'folded' | 'lognormal'
%   node_count  节点数（每个节点一条独立噪声）
%   steps       需要的样本数（= ceil(T/dt)+1，由 simulate_network 给）
%
% 输出
%   noise       node_count × steps，单位 pulses/s
%
% 做什么
%   把"不同阶段用不同均值噪声"这件实验性的事收成一个参数：
%     - 只给一个均值 → 等价于原来的"全程同分布高斯噪声"；
%     - 给 noise_seg_means = [90 60 90] → 三段均值不同的拼接噪声；
%     - 改 noise_segments = [0.5 0.5] → 两段，比例任意（会自动归一化）。
%   分段拼接本身没有生理含义，只是为了让一次长仿真里 p(t) 的统计量随时间变化，
%   用于观察增益动力学在"驱动强度切换"下的响应。论文正文不应把它写成模型的一部分。
%
% 与旧实现的数值一致性
%   旧版在 simulate_network 里写死三段 [0.3 0.5 0.2]，逐段调用 gauss() 再拼接。
%   本函数默认参数下逐段调用的次数、时长、顺序完全相同，因此随机数消耗序列也
%   相同，输出与旧版逐位一致（已用 fig1/fig2/fig3 校验和比对，见 log）。
%
% 调用关系
%   上游：simulate_network（noise_type = 'gauss'）
%   下游：gauss（单个分段的分布抽样）
T   = cfg.simulation_time;
dt  = cfg.dt;

frac   = getfield_default(cfg, 'noise_segments', [0.3 0.5 0.2]);
frac   = frac(:)' / sum(frac);
ns     = numel(frac);
mu     = expand_per_segment(getfield_default(cfg, 'noise_seg_means',  []), ...
                            getfield_default(cfg, 'noise_mean',  90),  ns, 'noise_seg_means');
sigma  = expand_per_segment(getfield_default(cfg, 'noise_seg_sigmas', []), ...
                            getfield_default(cfg, 'noise_sigma', 30),  ns, 'noise_seg_sigmas');
method = getfield_default(cfg, 'noise_method', 'normal');

noise = zeros(node_count, steps);

for n = 1:node_count
    seg = cell(1, ns);
    for k = 1:ns
        % 每段多取最后一个样本，拼接时丢掉，保证段与段之间不留重叠样本
        seg{k} = gauss(T * frac(k), dt, mu(k), sigma(k), method);
    end

    % 拼接：前 ns-1 段各丢最后一个样本，最后一段按需截取
    parts = cell(1, ns);
    for k = 1:ns-1
        parts{k} = seg{k}(1:end-1);
    end
    parts{ns} = seg{ns};

    sig = [parts{:}];
    if numel(sig) >= steps
        noise(n, :) = sig(1:steps);
    else
        % 段比例之和不足 1 或 dt 不能整除时样本会偏少：按最后一段的分布续抽，
        % 旧实现这里是用 0 补齐（会让 p(t) 突然掉到 0），改成续抽更稳。
        extra = gauss((steps - numel(sig)) * dt, dt, mu(end), sigma(end), method);
        noise(n, :) = [sig, extra(1:(steps - numel(sig)))];
    end
end

if size(noise, 2) ~= steps
    error('噪声长度 %d 与仿真步数 %d 不一致', size(noise, 2), steps);
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

% ======================================================================
function v = expand_per_segment(v, fallback, ns, name)
% 空 → 所有段用 fallback；标量 → 广播；向量 → 长度必须等于段数
if isempty(v)
    v = repmat(fallback, 1, ns);
elseif isscalar(v)
    v = repmat(v, 1, ns);
elseif numel(v) == ns
    v = v(:)';
else
    error('%s 长度 %d 与分段数 %d 不一致（cfg.noise_segments 有 %d 段）', ...
          name, numel(v), ns, ns);
end
end
