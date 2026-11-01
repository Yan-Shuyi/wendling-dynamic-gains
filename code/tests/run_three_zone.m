function run_three_zone()
% RUN_THREE_ZONE  三区网络（SOZ / PZ / NIZ）第一次自发启动：看跨区次序。
%
%   run_three_zone()
%
% 输入
%   无（参数写在下面的 cases 表里）
%
% 输出
%   tests/results/three_zone.txt   逐例汇总表
%   tests/results/three_zone.mat   逐例 EEG / 增益轨迹 / 标签（供 plot_three_zone 画图）
%
% 设计
%   启动器：SOZ 的 gain_bias = 0.005（唯一一个有背景放电基线的节点）。理由见
%   log/2026-11-01_三区初始值与周期性.md：静息点 x 恒为 0，没有 bias 就永远不启动；
%   bias 与节点参数无关，谁有 bias 谁的 Ge 都爬到 4.25，所以只能给致痫区。
%   A 一律 3.25：A 才是"有没有慢节律"的开关，SOZ 若一上来就 4.60，"起始未发生"
%   就不成立。B/G 按 (B,G) 平面落位（A=4.60 那张图的 ICTAL 中心 / preONSET / 深处）。
%   耦合：C = 15（两节点驱动—响应标定），选择性交给 W 而不是 C。
%
% 四个对照
%   1 SOZ_only    单节点，只有启动器     —— 自身时程（无任何耦合输入）
%   2 net         三区耦合               —— 主结果：谁先、滞后多少
%   3 net_noNIZ   SOZ→NIZ 权重置 0       —— NIZ 是否只能经 SOZ 进来
%   4 net_noBias  三区结构但无启动器     —— 阴性对照：没有 SOZ 就谁都不动
%
% 备注
%   不用 ESF（esf_enable = false）：事件必须由慢动力学自己长出来，
%   否则"起始"就是写进输入的时刻，不是结论。

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); setup_paths();
resdir = fullfile(root, 'tests', 'results');

T = 900;            % 单例时长 (s)
C = 15;             % 全局耦合常数（calib_coupling_c 定）

% 区：      A     B    G    bias
SOZ = [3.25   22    8   0.005];
PZ  = [3.25   32   12   0     ];
NIZ = [3.25   50   22   0     ];

% W(i,j) = j → i；延迟 D(i,j)，单位 s
W_full  = [0 0 0; 1 0 0; 0.15 0 0];
W_noNIZ = [0 0 0; 1 0 0; 0    0 0];
D_full  = [0 0 0; 0.015 0 0; 0.025 0 0];

cases = { ...
    'SOZ_only',  [SOZ],                     [],       [],      1; ...
    'net',       [SOZ; PZ; NIZ],             W_full,   D_full,  3; ...
    'net_noNIZ', [SOZ; PZ; NIZ],             W_noNIZ,  D_full,  3; ...
    'net_noBias',[SOZ(1:3) 0; PZ; NIZ],      W_full,   D_full,  3};

names = {'SOZ','PZ','NIZ'};
fid = fopen(fullfile(resdir, 'three_zone.txt'), 'w');
fprintf(fid, '# 三区网络自发启动，T=%d s，C=%d，legacy + reset_tol=0.05，无 ESF\n', T, C);
fprintf(fid, '# case        zone  rms_base  rms_max  onset   Ge_max  Ge_tmax  Ge_end  Gs_end  Gf_end  labelmix\n');
fprintf(fid, '# rms 均为 0.5 Hz 高通后的 1 s 窗，与判据口径一致\n');

S = struct();
for k = 1:size(cases,1)
    tag = cases{k,1}; P = cases{k,2}; W = cases{k,3}; D = cases{k,4};
    n = cases{k,5};
    cfg = default_config();
    cfg.node_count = n;
    cfg.A = P(:,1)'; cfg.B = P(:,2)'; cfg.G = P(:,3)';
    cfg.gain_bias = P(:,4)';
    cfg.model_dynamic = 1;
    cfg.gain_law = 'legacy';
    cfg.gain_reset_tol = 0.05;
    cfg.esf_enable = false;
    cfg.simulation_time = T;
    cfg.resample_dy = true;
    cfg.seed = 20261101;
    cfg.coupling_gain = C;
    cfg.sc_matrix = W;
    cfg.delay_matrix = D;

    [y, dy, cfg] = simulate_network(cfg);
    fs = cfg.eeg_fs;
    t = (0:size(y,2)-1) / fs;

    % ---- 1 s 窗 rms（先 0.5 Hz 高通去直流）----
    % 放电期 EEG 有一个 20 mV 量级的去极化平台，直接算 rms 得到的是平台高度，
    % 不是振荡幅度；判据（classify_params_v2.hp_cut）也是先去直流，这里必须对齐。
    yh = zeros(size(y));
    for i = 1:n
        yh(i,:) = highpass(y(i,:)', 0.5, fs)';
    end
    w = round(1*fs); nr = floor(size(yh,2)/w);
    rr = zeros(n, nr);
    for i = 1:nr
        seg = yh(:,(i-1)*w+(1:w));
        rr(:,i) = sqrt(mean(seg.^2, 2));
    end
    tr = (0:nr-1) + 0.5;

    % ---- onset：首次持续 5 s 越过 max(3×基线, 1.0 mV) ----
    base = median(rr(:, tr <= 30), 2);
    onset = nan(n,1);
    for i = 1:n
        above = rr(i,:) > max(3*base(i), 1.0);
        q = 0;
        for j = 1:nr-5
            if all(above(j:j+4)), q = j; break; end
        end
        if q > 0, onset(i) = tr(q) - 0.5; end
    end

    % ---- v2 判据逐窗标签 ----
    p = classify_params_v2();
    lab_txt = strings(1, n);
    onsets_frac = zeros(1, n);
    tcell = cell(1, n); lcell = cell(1, n);
    for i = 1:n
        ref = make_reference(y(i, 1:min(end, round(10*fs))), fs, p);
        [tc, labs] = classify_timeseries_v2(y(i,:), fs, p, ref);
        tcell{i} = tc; lcell{i} = labs;
        nw = max(1, numel(labs));
        cnt = sum(strcmp(labs, 'ONSET'));
        onsets_frac(i) = cnt / nw;
        lab_txt(i) = sprintf('I%.2f/p%.2f/O%.2f/u%.2f', ...
            sum(strcmp(labs,'ICTAL'))/nw, sum(strcmp(labs,'preONSET'))/nw, ...
            cnt/nw, sum(strcmp(labs,'uncertain'))/nw);
    end

    % ---- 汇总 ----
    fprintf('\n=== %s (n=%d) ===\n', tag, n);
    fprintf('%-6s %8s %8s %8s %8s %8s %8s   %s\n', 'zone','rms_b','rms_mx','onset','Ge_mx','Ge_end','Ge_tmx','label mix');
    for i = 1:n
        Ge = dy(3*(i-1)+1, :);
        Gs = dy(3*(i-1)+2, :);
        Gf = dy(3*(i-1)+3, :);
        [gem, gi] = max(Ge);
        fprintf('%-6s %8.2f %8.2f %8.1f %8.2f %8.2f %8.1f   %s\n', ...
            names{i}, base(i), max(rr(i,:)), onset(i), gem, Ge(end), t(gi), char(lab_txt(i)));
        fprintf(fid, '%-12s %-5s %9.3f %8.2f %7.1f %7.2f %8.1f %7.2f %7.2f %7.2f  %s\n', ...
            tag, names{i}, base(i), max(rr(i,:)), onset(i), gem, t(gi), ...
            Ge(end), Gs(end), Gf(end), lab_txt(i));
    end
    % 注意：上面 onset 列用的是"3×基线 + 1 s 窗"，三区基线差 3 倍（1.00 / 0.31 / 0.38），
    % 同一相对阈值在不同区等价于不同绝对门槛，会得出"PZ 比 SOZ 早"的假象。
    % 跨区次序一律以 analyze_three_zone 的两层判据（t_elec / t_Ge）为准。
    if n > 1
        fprintf('  （跨区次序见 analyze_three_zone 的 t_elec / t_Ge，勿用上表 onset 列）\n');
    end

    e.y = y; e.yh = yh; e.dy = dy; e.t = t; e.rr = rr; e.tr = tr;
    e.base = base; e.onset = onset; e.labs = lcell; e.lt = tcell;
    e.frac = onsets_frac; e.names = names(1:n); e.T = T; e.C = C;
    S.(tag) = e;
end
fclose(fid);
save(fullfile(resdir, 'three_zone.mat'), '-struct', 'S');
fprintf('\n已写 results/three_zone.txt / .mat\n');
end
