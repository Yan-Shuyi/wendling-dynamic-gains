function check_noise_time_invariance()
% CHECK_NOISE_TIME_INVARIANCE  验证"仿真时长不该改变早期轨迹"。
%
%   check_noise_time_invariance()
%
% 为什么
%   旧默认 noise_segments = [0.3 0.5 0.2] 是按"占总时长的比例"切段的：T 从 600
%   改成 900，段边界从 180/480 挪到 270/720，同一个种子抽出来的噪声就换了一套，
%   于是"同一初始条件、同一参数"的早期轨迹也跟着变。这不是数值误差，是实验设计
%   被时长污染。
%
% 判据
%   单段噪声（默认）下：T = 300 与 T = 900 的前 300 s 应当逐点相同（差为 0）。
%   三段噪声下：不应当相同——用来确认上面那条不是"因为都没跑出东西"。
%
% 输出
%   tests/results/noise_time_invariance.txt
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); setup_paths();
resdir = fullfile(root, 'tests', 'results');

cfg0 = default_config();
cfg0.node_count = 1;
cfg0.A = 3.25; cfg0.B = 22; cfg0.G = 10;
cfg0.gain_bias = 0.005;
cfg0.model_dynamic = 1; cfg0.gain_law = 'legacy'; cfg0.gain_reset_tol = 0.05;
cfg0.esf_enable = false;
cfg0.seed = 20261101;

fid = fopen(fullfile(resdir, 'noise_time_invariance.txt'), 'w');
fprintf(fid, '# 同一 seed、同一参数，只改 T；比较前 300 s 的逐点最大差\n');
fprintf(fid, '# case              max|dEEG| (mV)   max|dGe| (mV)   identical\n');

for variant = {'single', 'three'}
    seg = 1; if strcmp(variant{1}, 'three'), seg = [0.3 0.5 0.2]; end
    Y = cell(1,2); G = cell(1,2); Ts = [300 900];
    for k = 1:2
        cfg = cfg0;
        cfg.simulation_time = Ts(k);
        cfg.noise_segments = seg;
        cfg.resample_dy = true;
        [y, dy] = simulate_network(cfg);
        n = round(300*cfg.eeg_fs);
        Y{k} = y(1, 1:n);
        G{k} = dy(1, 1:round(300*cfg.dy_fs));
    end
    dy_eeg = max(abs(Y{1} - Y{2}));
    dy_g   = max(abs(G{1} - G{2}));
    fprintf('%-18s %14.3e %15.3e   %s\n', variant{1}, dy_eeg, dy_g, ...
        string(dy_eeg == 0 && dy_g == 0));
    fprintf(fid, '%-18s %14.3e %15.3e   %s\n', variant{1}, dy_eeg, dy_g, ...
        string(dy_eeg == 0 && dy_g == 0));
end
fclose(fid);
fprintf('\n已写 results/noise_time_invariance.txt\n');
end
