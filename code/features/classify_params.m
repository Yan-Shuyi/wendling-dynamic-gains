% ==============================
% plaintext: 参数配置逻辑
% 1. 仿真基础参数：对齐你的仿真设置，Fs=1000，6s总时长，弃前1s瞬态
% 2. 扫描参数：固定A，扫描B(0~50mV)、G(0~30mV)
% 3. 一级分类阈值（幅值粗筛）：
% 核心修改：用「均值+极差+最小幅值」三重条件锁定低幅组
% 新增阈值说明：
%   - thresh_mean_low: 均值绝对值阈值，<此值才可能进低幅组（沿用之前的2mV经验值）
%   - thresh_min_allowed: 允许的最小幅值，>此值才是“正常下冲”，<此值说明有强棘波
%   - thresh_range_small: 小极差阈值，<此值说明幅值波动极小，判定为类型1
% 4. 低幅组[1,2,3]专用阈值：
%    - thresh_diff_ratio: 棘波候选占比阈值，低于此值为Type1（无棘波）
%    - thresh_peak_num: 频谱显著峰数量阈值，≤2为Type3，>2为Type2
% 5. 高幅组[4,5,6]专用阈值：
%    - thresh_high_freq_ratio: 15Hz以上能量占比阈值，高于此值为Type5
%    - thresh_spec_bw: 频谱90%能量带宽阈值，低于此值为Type6（纯净正弦），高于为Type4
% ==============================
function params = classify_params()
% CLASSIFY_PARAMS  Thresholds and window settings for the EEG state classifier.
%
%   params = classify_params()
%
% 输入
%   无
%
% 输出
%   params   结构体，含两级判定用到的全部阈值与滑窗设置。关键字段：
%              win_len / win_step     滑窗长度与步长（s），论文用 5 s / 1 s
%              drop_head              丢弃开头多少秒的瞬态
%              thresh_mean            一级判定的均值阈值（mV）：>= 此值进高幅组
%              thresh_range_low       低幅组内"极差小"的阈值（mV）
%              thresh_diff_max        低幅组内"最大一阶差分"阈值（mV）
%              thresh_peak_num        显著谱峰数阈值
%              thresh_hf_ratio        高频（>=15 Hz）能量占比阈值
%              thresh_bw90            90% 能量带宽阈值（Hz）
%
% 做什么
%   集中存放 §3.3 决策树的全部阈值，训练/校准与判定共用一份，避免两处写死
%   导致图与正文对不上。
%
% 调用关系
%   上游：run_fig2.m、classify_eeg_timeseries、classify_single
%   下游：无
%
% 备注
%   论文 Fig. 2 那棵树的每个判据数字都直接取自本文件；改阈值必须同步改论文。
%   文件内的注释写着"极差 2 mV、差分 3*std"，但实际生效的是 3.0 与 1.5，
%   以代码为准（论文已按代码写）。这处注释与代码的出入见 log。
params.Fs = 256;                  % 采样率
params.total_duration = 6;          % 单组参数总仿真时长(s)
params.transient_duration = 1;       % 舍弃前1s瞬态数据(s)
params.noise_mu = 90;                % 输入高斯噪声均值
params.noise_var = 30;               % 输入高斯噪声方差
params.A_fixed = 5;                   % 固定A值(mV)，可修改为3~7mV扫描不同A的热力图

params.win_half = 2.5;             % 滑窗半宽(秒)，±2.5s共5s分析窗
params.step_sec = 1;               % 滑窗步长(秒)，1秒一步兼顾效率和分辨率

% 参数扫描范围（完全对齐第二篇文献）
% ----------------------
params.B_list = 0:1:50;             % 慢抑制增益B扫描范围(mV)，步长1mV
params.G_list = 0:1:30;             % 快抑制增益G扫描范围(mV)，步长1mV

%% ----------------------
% 2. 【核心】一级分类阈值（均值一刀切，无交叉）
% ----------------------
params.thresh_mean_high = 5.0;     % 均值阈值(mV)：>=5mV→高幅组[4,5,6]；<5mV→低幅组[1,2,3]
                                   % 校准建议：跑典型Type4的最小均值，加0.5mV缓冲

%% ----------------------
% 3. 低幅组[1,2,3]专用阈值（仅在mean<5mV时生效，负值仅在此处讨论）
% ----------------------
params.thresh_range_123 = 3.0;     % 极差阈值(mV)：<2mV→Type1（平稳背景）；>=2mV→进入棘波筛查
params.diff_max = 1.5;          % 差分阈值：diff_max > 2→判定为棘波陡变
params.thresh_peak_num = 2;         % 频谱尖峰数阈值：<=2→Type3（持续棘波，频谱纯净）；>2→Type2（零星棘波，背景杂乱）
params.min_peak_height_ratio = 0.1;% 频谱显著峰最小高度：≥最大功率的10%（过滤噪声假峰）
params.min_peak_dist = 2;          % 频谱峰最小间隔(Hz)：≥2Hz（避免相邻噪声峰被误判）

%% ----------------------
% 4. 高幅组[4,5,6]专用阈值（仅在mean>=5mV时生效，绝不检查最小值/负值）
% ----------------------
params.high_freq_start = 15;       % 高频起始频率(Hz)：15Hz以上算高频能量
params.thresh_high_freq_ratio = 0.3;% 高频能量占比阈值：>30%→Type5（低幅快活动）
params.thresh_spec_bw = 2.0;       % 频谱90%能量带宽阈值(Hz)：<2Hz→Type6（高幅准正弦，带宽极窄）；>=2Hz→Type4（慢节律，带宽较宽）

%% ----------------------
% 5. 特征维度（无需修改，对应extract_features的输出顺序）
% ----------------------
params.feature_dim = 11; % [mean_val, range_val, diff_max, f_peak, high_freq_ratio, peak_num, spec_bw, min_val, max_val, rms, line_length]

%% ----------------------
% 6. 可视化配置
% ----------------------
params.interp_to_signal = true;   % 是否将类型序列插值到和原始信号等长（方便对齐）
params.type_colors = {             % 6类配色（和热力图完全统一）
    [0 0 1],   % 1: 蓝（正常背景）
    [0 1 1],   % 2: 青（零星棘波）
    [0 1 0],   % 3: 绿（持续棘波）
    [1 0.5 0], % 4: 橙（慢节律）
    [1 0 0],   % 5: 红（低幅快活动）
    [0.85 0.85 0.85] % 6: 浅灰（高幅准正弦，避免白色看不见）
};
params.type_names = {'1:Normal','2:Sporadic','3:Sustained','4:Slow','5:Fast','6:High-amp'};

end