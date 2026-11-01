function [outputs_fs, dy_fs, cfg] = simulate_network(cfg)
% SIMULATE_NETWORK  Run one simulation of the (possibly decoupled) Wendling network.
%
%   [outputs_fs, dy_fs, cfg] = simulate_network(cfg)
%
% 输入
%   cfg   由 default_config() 得到，调用方按需覆盖字段
%
% 输出
%   outputs_fs  重采样后的 EEG，node_count × n，单位 mV
%   dy_fs       重采样后的 E/S/F 轨迹，3*node_count × n；
%               cfg.resample_dy = false 时返回 []
%   cfg         回填后的配置（补全 steps、实际的 A/B/G 向量等，便于复现）
%
% 做什么
%   1. 按 cfg 逐节点构造 wendling 模型对象（A/B/G 若是向量则逐节点取值）；
%   2. 组装结构连接矩阵与延迟矩阵（当前论文三张图都是全零 = 完全解耦）；
%   3. 构造 System 对象，逐步推进网络动力学，记录 EEG 与增益轨迹；
%   4. 生成外部驱动 p(t)：默认分段高斯噪声（noise_type = 'gauss'，见
%      noise/gauss_noise.m）；实验分支用另一个 Wendling 输出当驱动
%      （noise_type = 'wendling'，见 noise/wendling_input.m）；
%   5. 把原始采样率重采样到 EEG / dy 的目标采样率。
%
% 调用关系
%   上游：run_fig1.m / run_fig2.m / run_fig3.m
%   下游：wendling（模型）、System（网络推进）、gauss_noise / wendling_input
%         （驱动源）、resample_fs（重采样）
%
% 备注
%   本函数由原 main.m / main_fig1.m / main_fig2.m / main_fig3.m 的仿真循环原样
%   抽取而来，循环体逐行保持一致，只把写死的差异改成读 cfg。
%   重构前后的 fig2 已做逐像素比对，输出完全一致（见 log）。
if ~isempty(cfg.seed)
    rng(cfg.seed);      % 固定随机种子，保证图可复现
end

N = cfg.node_count;

%% ---------------- 参数展开为逐节点向量 ----------------
A = expand_param(cfg.A, N);
B = expand_param(cfg.B, N);
G = expand_param(cfg.G, N);
cfg.A = A; cfg.B = B; cfg.G = G;

%% ---------------- 模型初始化 ----------------
models = cell(1, N);
for i = 1:N
    models{i} = wendling(A(i), B(i), G(i), cfg.dt, cfg.model_dynamic);
    % 增益复位容差（长时程"发作—恢复—再发作"需要；默认 0 = 旧行为）
    if isfield(cfg, 'gain_reset_tol') && ~isempty(cfg.gain_reset_tol)
        models{i}.gain_reset_tol = cfg.gain_reset_tol;
    end
    % 背景放电基线（可逐节点给；默认 0 = 旧行为）
    if isfield(cfg, 'gain_bias') && ~isempty(cfg.gain_bias)
        gb = cfg.gain_bias;
        if isscalar(gb), gb = repmat(gb, 1, N); end
        models{i}.gain_bias = gb(i);
    end
    % 增益动力学形式（D13）：'legacy' 为默认，'bistable' 见 model/gain_bistable.m
    if isfield(cfg, 'gain_law') && ~isempty(cfg.gain_law) ...
            && ~strcmpi(cfg.gain_law, 'legacy')
        models{i} = models{i}.set_gain_law(cfg.gain_law, cfg.bs_params);
    end
end

%% ---------------- 系统维度 ----------------
cfg.state_size  = sum(cellfun(@(m) m.state_size,  models));
cfg.dy_size     = sum(cellfun(@(m) m.dy_size,     models));
cfg.input_size  = sum(cellfun(@(m) m.input_size,  models));
cfg.output_size = sum(cellfun(@(m) m.output_size, models));

%% ---------------- 网络与噪声 ----------------
if isempty(cfg.sc_matrix)
    sc_matrix = zeros(N, N);
else
    sc_matrix = cfg.sc_matrix;
end
sc_matrix(eye(N) == 1) = 0;                 % 确保不自连

if isempty(cfg.delay_matrix)
    delay_matrix = zeros(N, N);
else
    delay_matrix = cfg.delay_matrix;
end

system  = System(models, cfg, sc_matrix, delay_matrix);
steps   = ceil(cfg.simulation_time / cfg.dt) + 1;
noise   = build_noise(cfg, steps);

%% ---------------- 可选的瞬时扰动（双稳增益律的触发方式） ----------------
% 旧版触发 = 直接把增益跳到 esf_params（阶跃）。双稳形式下更自然的是给一段
% 短暂的外部驱动：活动越过阈值后由双稳项自己把增益闩住，之后即使扰动消失
% 增益也不会回落（滞后）。两者用 cfg.esf_gain_jump / cfg.esf_pulse 分别开关。
if cfg.esf_enable && isfield(cfg, 'esf_pulse') && ~isempty(cfg.esf_pulse) ...
        && cfg.esf_pulse ~= 0
    i1 = round(cfg.esf_time / cfg.dt) + 1;
    i2 = min(steps, i1 + round(cfg.esf_pulse_duration / cfg.dt) - 1);
    noise(:, i1:i2) = noise(:, i1:i2) + cfg.esf_pulse;
end

%% ---------------- 仿真主循环 ----------------
inputs   = zeros(cfg.input_size,  steps);
dy_param = zeros(cfg.dy_size,     steps);
outputs  = zeros(cfg.output_size, steps);

esf_step = round(cfg.esf_time / cfg.dt);
gain_jump = true;
if isfield(cfg, 'esf_gain_jump') && ~isempty(cfg.esf_gain_jump)
    gain_jump = logical(cfg.esf_gain_jump);
end

for i = 1:steps
    [~, dy_param(:, i), inputs(:, i), outputs(:, i)] = system.network_dynamics(noise(:, i));

    if cfg.esf_enable && gain_jump && i == esf_step
        for j = 1:N
            if cfg.esf_jitter > 0
                ge = cfg.esf_params(1) + cfg.esf_jitter * randn;
            else
                ge = cfg.esf_params(1);
            end
            system.models{j}.ESF(ge, cfg.esf_params(2), cfg.esf_params(3));
        end
    end
end

%% ---------------- 重采样 ----------------
outputs_fs = resample_fs(outputs, 1/cfg.dt, cfg.eeg_fs, cfg.use_parallel);

if cfg.resample_dy && cfg.dy_size > 0
    dy_fs = resample_fs(dy_param, 1/cfg.dt, cfg.dy_fs, cfg.use_parallel);
else
    dy_fs = [];      % 静态模式 dy_size = 0，没有轨迹可重采样
end
end

% ======================================================================
function p = expand_param(v, N)
% 标量参数复制到所有节点；长度不符的向量直接报错，避免"静默用错参数"
if isscalar(v)
    p = repmat(v, 1, N);
elseif numel(v) == N
    p = v(:)';
else
    error('参数长度 %d 与节点数 %d 不一致', numel(v), N);
end
end

% ======================================================================
function noise = build_noise(cfg, steps)
% 输入噪声分发：把 cfg.noise_type 翻译成一个具体的驱动源。
%   'gauss'    分段高斯白噪声（论文三张图都用这个）
%   'wendling' 另一个 Wendling 模型的输出（实验分支，见 noise/wendling_input.m）
% 兼容旧的数字写法：0 → 'gauss'，1 → 'wendling'（旧版 1 是直接报错的占位分支）。
% 分段方式（每段多长、各段均值/标准差）现在全部是参数，见 gauss_noise.m 头注释。
N = cfg.node_count;

nt = cfg.noise_type;
if isnumeric(nt)
    switch nt
        case 0, nt = 'gauss';
        case 1, nt = 'wendling';
        otherwise, error('未知 noise_type = %d', nt);
    end
end

switch lower(nt)
    case 'gauss'
        noise = gauss_noise(cfg, N, steps);
    case 'wendling'
        noise = wendling_input(cfg, N, steps);
    otherwise
        error('未知 noise_type = ''%s''（可选 ''gauss'' / ''wendling''）', nt);
end

if size(noise, 2) ~= steps
    error('噪声长度 %d 与仿真步数 %d 不一致', size(noise, 2), steps);
end
end
