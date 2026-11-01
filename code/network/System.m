classdef System < handle
% SYSTEM  把 N 个模型对象推进一个时间步的网络层（含延迟与结构耦合）
%
%   sys = System(models, config, sc_matrix, delay_matrix)
%   sys = System(models, config, sc_matrix, delay_matrix, sig)
%   [x, dy_param, inputs, outputs] = sys.network_dynamics(noise)
%
% 输入
%   models         cell 数组，每个元素是一个 base 子类对象（wendling）
%   config         default_config() 的结构体（用到 node_count / dt /
%                  state_size / dy_size / output_size）
%   sc_matrix      N×N 结构连接矩阵；第 (i,j) 元 = 从 j 到 i 的耦合强度
%                  （全零 = 节点完全解耦，当前论文三张图都是全零）
%   delay_matrix   N×N 传导延迟矩阵（秒）；第 (i,j) 元 = 信号从 j 传到 i
%                  花费的时间，与 sc_matrix(i,j) 一一对应
%   sig            （可选）S 型函数的三个参数，结构体字段 e0 / r / v0
%                  省略或给 [] 时自动从 models{1} 读取同名属性
%
% 输出
%   x         拼接后的全部节点状态向量（N × state_size）
%   dy_param  拼接后的增益导数（N × dy_size）
%   inputs    本步每个节点实际收到的输入（噪声 + 耦合项）
%   outputs   本步每个节点的模型输出（即 EEG 采样值）
%
% 做什么
%   1. 从环形延迟缓冲区里按 delay_matrix(i,j) 取回"若干步之前"的节点 j 输出；
%   2. 先过 S 型函数、再用 sc_matrix(i,j) 加权累加到节点 i 的输入；
%   3. 逐个节点调用 model.simulate()，把状态/导数/输出拼成整块向量；
%   4. 把本步输出写回延迟缓冲区，指针环形前进。
%
% 耦合顺序（论文里必须和这里一致）
%   生效的是"先 S 再耦合"，并带一个全局耦合常数 C：
%       p_i(t) = noise_i + C · Σ_j SC(i,j) · S(y_out,j(t - D(i,j)))
%   被注释掉的另一种写法是"先耦合再 S"：
%       p_i(t) = noise_i + C · S( Σ_j SC(i,j) · y_out,j(t - D(i,j)) )
%   两者不等价。论文正文应写成 C · Σ_j W_ij · S(y_out,j)，与代码一致。
%   C 由 config.coupling_gain 给，默认 1（= 不乘 C，存量图数值不变）。
%
% 负向是否参与耦合（config.coupling_clip_negative）
%   旧实现把 y_out < 0 一律钳成 0 再进 S，等于"抑制性 PSP 不参与耦合"。
%   后果：在噪声均值较高、y_out 长期为负的参数下，耦合项退化成恒定的
%   S(0) = 2*e0/(1+exp(r·v0)) ≈ 0.168，D 怎么改都看不到变化（"耦合假死"）。
%   现在把它做成开关：
%     true （默认，向后兼容）—— 钳到 0，与已发表的三张图完全一致；
%     false —— 负值照常进 S，S 对负输入单调、取值很小但不为零，
%              抑制性投射因此能对下游产生"减法"作用，生理上更合理。
%   注意：设成 false 会改变 fig2 / fig3 的数值，切换后必须重跑并重画。
%
% 已知隐患（论文结束前不必动，但换参数/换网络前必须正面处理）
%   1. D1：S 的饱和值 2*e0 = 5，而噪声均值 noise_mean = 90，两者相差近 20 倍。
%      耦合项在 SC 取 DTI 量级（<0.1）时只有 0.5 左右，等于被噪声完全淹没。
%      解法已落到 coupling_gain：标定 C 使 C·max(ΣW)·2*e0 与噪声同量级。
%      量级参考 C ≈ noise_mean / (2*e0·max列和 W)（噪声 90、列和≈1 时 C≈18）。
%   2. 延迟缓冲一律从 1 开始读（缓冲未填满时用最早的一格），等价于"开机
%      瞬间所有节点初值相同"，因此前 max_delay_steps 步的耦合是伪信号。
%      单节点或零延迟时无影响；做长延迟网络时应丢弃前 max(D(:)) 秒。
%   3. 每条连接的延迟只取 D(i,j) 的一个值，不支持分布式延迟（一张 SC 一个 D）。
%
% 调用关系
%   上游：simulate_network（每步调用一次）
%   下游：models{i}.simulate
    properties
        models          % 存储 N 个节点的模型对象
        config          % 仿真配置（dt、节点数、维度等）
        sc_matrix       % 结构连接矩阵
        delay_matrix    % 延迟矩阵（秒）
        delay_history   % 延迟历史矩阵：节点数 × 最大延迟步数
        delay_ptr       % 每个节点的写入指针（位置）
        max_delay_steps % 最大延迟步数（由最大延迟秒数/dt得到）
        sig             % S 型函数参数：e0 / r / v0
        S_fun           % 已绑定参数的 S 型函数句柄（避免每步重建）
        coupling_gain   % 全局耦合常数 C（config.coupling_gain，默认 1）
        clip_negative   % 耦合前是否把 y_out<0 钳成 0（config.coupling_clip_negative）
        results         % 结果存储
    end

    methods
        % 构造函数：初始化系统
        function obj = System(models, config, sc_matrix, delay_matrix, sig)
            obj.models = models;                  % 传入模型列表
            obj.config = config;                  % 传入配置
            obj.sc_matrix = sc_matrix;            % 结构连接矩阵
            obj.delay_matrix = delay_matrix;      % 延迟矩阵（秒）
            obj.results = struct();               % 结果结构体

            % ---- 全局耦合常数 C 与负向钳位开关（缺省 = 旧行为）----
            % coupling_gain 的作用：把结构矩阵的整体尺度与"节点间同步强度"解耦。
            % D1 问题（S 饱和于 2*e0 = 5，而噪声均值 90）就靠它标定：
            % 取 C ≈ noise_mean / (2*e0*max(ΣW)) 量级，才能让耦合项与噪声同量级。
            % 默认 1，等于不乘 C，三张存量图数值不变。
            obj.coupling_gain = get_cfg(config, 'coupling_gain', 1);
            obj.clip_negative = get_cfg(config, 'coupling_clip_negative', true);
            if ~isscalar(obj.coupling_gain) || ~isnumeric(obj.coupling_gain)
                error('coupling_gain 必须是标量，当前为 %s', mat2str(size(obj.coupling_gain)));
            end

            % ---- S 型函数参数：默认跟随模型，也允许调用方覆盖 ----
            if nargin < 5 || isempty(sig)
                m1 = models{1};
                if isprop(m1, 'e0') && isprop(m1, 'r') && isprop(m1, 'v0')
                    sig = struct('e0', m1.e0, 'r', m1.r, 'v0', m1.v0);
                else
                    % 模型没有这三个属性时的兜底值（等于 Wendling 原文标准值）
                    sig = struct('e0', 2.5, 'r', 0.56, 'v0', 6);
                end
            end
            obj.sig = sig;

            % 把数值先取成局部变量再建句柄：句柄一旦捕获 obj 本身，
            % handle 类会出现自引用，对象永远不会被释放。
            e0 = obj.sig.e0;  r = obj.sig.r;  v0 = obj.sig.v0;
            obj.S_fun = @(v) 2*e0 ./ (1 + exp(r*(v0 - v)));

            dt = config.dt;                      % 仿真步长
            max_delay_sec = max(delay_matrix(:)); % 所有延迟里最大的秒数
            obj.max_delay_steps = ceil(max_delay_sec / dt) + 10; % 转成步数，多给10格余量

            N = config.node_count;                % 节点数
            obj.delay_history = zeros(N, obj.max_delay_steps); % 初始化历史矩阵
            obj.delay_ptr = ones(N, 1);           % 每个节点从第1位开始写
        end

        % ================= 核心：全矩阵化网络动力学函数 ================
        function [x, dy_param, inputs, outputs] = network_dynamics(obj, noise)
            N     = obj.config.node_count;         % 节点数
            dt    = obj.config.dt;                 % 步长
            D     = obj.delay_matrix;              % 延迟矩阵（秒）
            SC    = obj.sc_matrix;                 % 结构连接矩阵
            M     = obj.max_delay_steps;           % 最大延迟步数
            ptr   = obj.delay_ptr;                 % 写入指针
            hist  = obj.delay_history;             % 延迟历史矩阵：N × M
            S     = obj.S_fun;                     % S 型函数（参数已在构造时绑定）

            % 初始化输出
            x       = zeros(obj.config.state_size, 1);
            dy_param = zeros(obj.config.dy_size, 1);
            outputs = zeros(obj.config.output_size, 1);

            %% ==================== 矩阵耦合（无自环） ====================
            SC_no_self = SC;
            SC_no_self(1:N+1:end) = 0;  % 对角线清0

            % ==================== 按 (i,j) 配对累加耦合输入 ====================
            % 对每一对"发送 j → 接收 i"，用 D(i,j) 的延迟取节点 j 的输出，
            % 再过 S 型函数、乘 SC(i,j)。旧版用 D(j,j)（矩阵对角线），
            % 只有当所有延迟相同时才碰巧正确；D 一旦非对称就会取错信号。
            % 解耦时（SC 全零）整个循环被跳过，input_vec 退化为纯噪声，
            % 与旧版完全等价，且不产生额外开销。
            input_vec = noise(:);
            if any(SC_no_self(:) ~= 0)
                C  = obj.coupling_gain;            % 全局耦合常数
                do_clip = obj.clip_negative;       % 是否把负向输出钳成 0
                for j = 1:N
                    for i = 1:N
                        w_ij = SC_no_self(i, j);
                        if w_ij == 0
                            continue;
                        end
                        steps_back = round(D(i, j) / dt);
                        pos = ptr(j) - 1 - steps_back;
                        if pos < 1
                            pos = 1;               % 缓冲未填满时用最早的一格
                        end
                        pos = min(pos, M);         % 绝不超过缓冲区长度

                        v = hist(j, pos);
                        if do_clip && v < 0
                            v = 0;                 % 旧行为：负向 PSP 不参与耦合
                        end
                        input_vec(i) = input_vec(i) + C * w_ij * S(v);
                    end
                end
            end

            %% ==================== 更新所有节点 ====================
            idx_s = 1;
            idx_dy = 1;
            idx_o = 1;
            y = zeros(N, 1);

            for i = 1:N
                m = obj.models{i};
                [dy_ESF, xi, yi] = m.simulate(input_vec(i));

                x(idx_s : idx_s + m.state_size - 1) = xi;
                dy_param(idx_dy : idx_dy + m.dy_size - 1) = dy_ESF;
                outputs(idx_o : idx_o + m.output_size - 1) = yi;
                y(i) = yi;

                idx_s = idx_s + m.state_size;
                idx_dy = idx_dy + m.dy_size;
                idx_o = idx_o + m.output_size;
            end
            inputs = input_vec;

            %% ==================== 写入当前值到延迟矩阵 ====================
            lin_idx = sub2ind([N, M], (1:N)', ptr);
            hist(lin_idx) = y;

            % 指针前进，环形
            ptr = ptr + 1;
            ptr(ptr > M) = 1;

            % 写回对象
            obj.delay_history = hist;
            obj.delay_ptr     = ptr;
        end

    end
end

% ======================================================================
function v = get_cfg(config, name, default)
% 取配置字段，缺省或为空时用默认值（保证老配置文件也能跑）
if isstruct(config) && isfield(config, name) && ~isempty(config.(name))
    v = config.(name);
else
    v = default;
end
end
