function cfg = default_config()
% DEFAULT_CONFIG  Build the default configuration struct for one simulation.
%
%   cfg = default_config()
%
% 输入
%   无
%
% 输出
%   cfg   结构体，字段含义见下方"字段表"；调用方按需覆盖个别字段即可
%
% 做什么
%   把一次仿真需要的全部参数收敛到一个结构体里。以前 main.m / main_fig1.m /
%   main_fig2.m / main_fig3.m 各自复制了约 200 行几乎相同的仿真代码，差异只有
%   节点数、是否动态增益、是否触发 ESF、是否重采样 dy 这四五项；现在这些差异
%   全部落成字段，仿真循环只保留一份（simulate_network.m）。
%
% 字段表（★ = 三个 run_fig*.m 里经常需要覆盖的）
%   ★node_count      仿真多少个互相独立的节点（各 run_fig*.m 自己设）
%   ★A, B, G         初始 (G_e, G_s, G_f)，单位 mV；B/G 可给向量，逐节点取值
%   ★model_dynamic   0 = 静态增益（经典 Wendling）；1 = 增益由 ESF 迟滞律驱动
%   ★esf_enable      是否在 esf_time 处触发一次增益跳变
%   ★esf_time        触发时刻（s）
%   ★esf_params      触发后跳到的 [G_e, G_s, G_f]
%   ★esf_jitter      触发增益的随机抖动幅度；0 = 每次都一样（单例图用）
%   ★resample_dy     是否需要把 E/S/F 轨迹重采样输出（只有画轨迹的图才要）
%    simulation_time  单次仿真时长（s），默认 200
%    dt               积分步长（s），1/256，与 Wendling 原文采样率一致
%    eeg_fs           EEG 输出的采样率（Hz），默认 256；dy_fs 为增益轨迹采样率
%    noise_type       驱动源类型：'gauss'（默认，分段高斯噪声）/ 'wendling'（实验分支）
%    noise_mean/sigma 外部传入驱动 p(t) 的高斯噪声均值/标准差（pulses/s）
%    noise_segments   分段比例，默认 1（单段）。分段按"占总时长的比例"切，
%                     所以 T 一变、段边界就变、同一个种子的早期轨迹也跟着变。
%                     单段时 T 只决定生成多少样本，前面的样本逐点不变。
%                     论文 Fig.1–5 的存量结果基于 [0.3 0.5 0.2]，重跑需显式指定。
%    noise_seg_means/sigmas  逐段均值/标准差；[] = 各段相同（旧版行为）
%    noise_wendler_*  noise_type = 'wendling' 时驱动模型的参数与量纲开关
%    sc_matrix        结构连接矩阵；全零 = 节点完全解耦（当前论文三张图都是）
%    delay_matrix     传导延迟矩阵（s）；本阶段搁置，保持 [] 即可
%    coupling_gain    全局耦合常数 C，实际权重 = C · sc_matrix；默认 1（等价于不乘 C）
%    coupling_clip_negative  耦合前是否把 y_out<0 钳成 0；true = 旧行为
%    use_parallel     重采样是否用 parfor
%
% 调用关系
%   下游：simulate_network
%   上游：run_fig1.m / run_fig2.m / run_fig3.m
%
% 备注
%   节点之间当前是"零耦合"：sc_matrix 全零，所以 node_count = N 实际等价于
%   "做 N 次独立实验"。论文里必须按独立重复来写，不能写成 N 节点耦合网络。
cfg = struct();

%% ---------------- 时间与采样 ----------------
cfg.simulation_time = 200;      % 仿真时长 (s)
cfg.dt              = 1/256;    % 积分步长 (s)
cfg.eeg_fs          = 256;      % EEG 重采样率 (Hz)
cfg.dy_fs           = 256;      % E/S/F 重采样率 (Hz)
cfg.resample_dy     = true;     % 是否需要 E/S/F 轨迹（fig1 静态扫描不需要）
cfg.use_parallel    = false;    % resample 是否并行

%% ---------------- 节点与参数 ----------------
cfg.node_count   = 1;           % 节点数
cfg.A            = 3.25;        % 兴奋性增益 (mV)；标量=所有节点同值，向量=逐节点
cfg.B            = 22;          % 慢抑制增益 (mV)
cfg.G            = 10;          % 快抑制增益 (mV)
cfg.model_dynamic = 1;          % 1 = 增益由 ESF 迟滞律驱动；0 = 增益固定（静态）

%% ---------------- 网络耦合 ----------------
cfg.sc_matrix    = [];          % [] 表示完全解耦（全零）
cfg.delay_matrix = [];          % [] 表示零延迟（延迟路径已实现，本阶段搁置不用）
cfg.coupling_gain = 1;          % 全局耦合常数 C：实际耦合权重 = C · sc_matrix(i,j)
                                % 见 System.m 头注释；取 1 时与"不乘 C"完全等价
                                % 耦合网络用 C = 15：由两节点驱动—响应扫描定（tests/calib_coupling_c.m），
                                % 此时一条满权重边贡献 16.7 pulses/s ≈ μ/5；C>20 时背景节点会被判据误判为发作
cfg.coupling_clip_negative = true;  % 耦合前是否把 y_out < 0 钳成 0
                                % true  = 旧行为（负向 PSP 不参与耦合）
                                % false = 负值照常进 S 函数（生理上更说得通，但会改图）

%% ---------------- 输入噪声 ----------------
cfg.noise_type   = 'gauss';     % 'gauss'  = 分段高斯白噪声（论文用这个）
                                % 'wendling' = 另一个 Wendling 模型的输出当驱动
                                % 兼容旧数字写法：0 → 'gauss'，1 → 'wendling'
cfg.noise_mean   = 90;          % 高斯噪声均值（pulses/s）
cfg.noise_sigma  = 30;          % 高斯噪声标准差
cfg.noise_segments   = 1;       % 各段占总时长的比例（会自动归一化）
                                % 1 = 单段（默认）：与 T 无关，早期轨迹不随时长变
                                % [0.3 0.5 0.2] = 论文 Fig.1–5 的历史设置，
                                %   它按比例切段 → T 改了就换了一套噪声，跨 T 不可比
cfg.noise_seg_means  = [];      % 逐段均值；[] = 全部用 noise_mean
                                % 例：给 [90 60 90] 即"三段不同均值"，此前只能改代码
cfg.noise_seg_sigmas = [];      % 逐段标准差；[] = 全部用 noise_sigma
cfg.noise_method     = 'normal';% 'normal' | 'truncated' | 'folded' | 'lognormal'

%% ---------------- 替代驱动源（noise_type = 'wendling' 时生效）----------------
cfg.noise_wendler_params    = [4.25 22 10]; % 驱动模型的 [G_e G_s G_f]
cfg.noise_wendler_dynamic   = 0;    % 驱动模型是否启用 ESF 增益动力学
cfg.noise_wendler_normalize = true; % 是否把驱动输出的 mV 映射到噪声量纲
                                    % false = 直接用原始 y_out（mV）驱动下游

%% ---------------- 病理触发（ESF） ----------------
cfg.esf_enable   = true;        % 是否在 t = esf_time 处改写突触增益
cfg.esf_time     = 10;          % 触发时刻 (s)
cfg.esf_params   = [4.25, 22, 10];  % 触发后 (Ge, Gs, Gf)
cfg.esf_jitter   = 0;           % Ge 的随机抖动幅度；fig3 用 0.05，fig2 用 0
cfg.esf_gain_jump  = true;      % 触发时是否直接把增益写到 esf_params
cfg.esf_pulse      = 0;         % 可选：触发后附加的外部驱动脉冲（pulses/s）
                                % 只给 bistable 增益律用（"瞬时扰动越过双稳阈值"）
cfg.esf_pulse_duration = 2;     % 该脉冲持续时长 (s)

%% ---------------- 增益动力学形式（D13） ----------------
cfg.gain_law  = 'legacy';       % 'legacy'   = 分段线性 + 三个状态机（存量三图用这个）
                                % 'bistable' = 统一三次型双稳方程（见 model/gain_bistable.m）
cfg.gain_reset_tol = 0;         % 增益复位容差（mV）。0 = 旧行为（要求精确回到静息值，
                                % 浮点下几乎不可能，于是 legacy 律下每个节点只发一次）；
                                % 取 0.05 才能在长时程里出现"恢复→再发作"。200 s 的存量图不受影响
cfg.gain_bias = 0;              % 背景放电基线（mV），加在增益驱动 x 上。
                                % 0 = 旧行为（静息 x 恒为 0，增益不动）；
                                % 取 0.002–0.005 时 G_e 会缓慢爬向 3.65–4.25 mV，
                                % 长时程里可自发进入发作并循环（三区网络用逐节点向量）
cfg.bs_params = [];             % 双稳参数覆盖；[] = 用 wendling 里的默认值
                                % 字段：Glo span phi0 phi_thr A B hc C0 C1，均为 3×1

%% ---------------- 复现性 ----------------
cfg.seed         = 20261029;    % [] 表示不固定随机种子

end
