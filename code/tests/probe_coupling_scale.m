function probe_coupling_scale()
% PROBE_COUPLING_SCALE  标定全局耦合常数 C 前的量级探测。
%
% 问的问题：耦合项是 C·W·S(y_out)。要定 C，先得知道 S(y_out) 在这个模型的
% 典型参数点上到底有多大——因为 S 的饱和值只有 2*e0 = 5（脉冲/秒），
% 而噪声均值是 90。两者差 18 倍，C 的量级完全由 S(y_out) 的实际分布决定。
%
% 做法：三个候选参数点（NIZ / PZ / SOZ）各跑 60 s 静态增益，
% 取 t>20 s，统计 y_out 与 S(y_out) 的均值、分位数、rms，并换算
% "C 使耦合项等于噪声均值"所需的 C*。
%
% 输出：命令行表格 + tests/results/calib_c_scale.txt

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); setup_paths();

pts = { ...
    'NIZ', 3.25, 45, 25; ...
    'PZ',  4.00, 35, 10; ...
    'SOZ', 4.60, 25, 10};

e0 = 2.5; r = 0.56; v0 = 6;
Sfun = @(v) 2*e0 ./ (1 + exp(r*(v0 - v)));

fid = fopen(fullfile(root, 'tests', 'results', 'calib_c_scale.txt'), 'w');
fprintf(fid, 'point  A     B   G    rms_y   y_med   y_p05  y_p95   S_mean  S_p50  S_p95  C*_eqnoise\n');

fprintf('\n%-5s %5s %5s %5s %8s %8s %8s %8s %8s %8s %8s %9s\n', ...
    'pt','A','B','G','rms','y_med','y_p95','S_mean','S_p50','S_p95','C*','');
for k = 1:size(pts, 1)
    cfg = default_config();
    cfg.node_count = 1;
    cfg.A = pts{k,2}; cfg.B = pts{k,3}; cfg.G = pts{k,4};
    cfg.model_dynamic = 0;          % 静态增益：量的是"这个工作点本身的 S 分布"
    cfg.esf_enable = false;
    cfg.simulation_time = 60;
    cfg.resample_dy = false;
    cfg.seed = 20261031;

    y = simulate_network(cfg);      % 1 × n (mV)
    y = y(:)';
    t = (0:numel(y)-1) / cfg.eeg_fs;
    y = y(t > 20);                  % 丢掉前 20 s

    sv = Sfun(y);
    ymed = median(y); yp05 = prctile(y,5); yp95 = prctile(y,95);
    Smean = mean(sv); Sp50 = median(sv); Sp95 = prctile(sv,95);
    rmsv = sqrt(mean(y.^2));
    Cs = cfg.noise_mean / Smean;    % 使 C·1·S_mean 恰等于噪声均值

    fprintf('%-5s %5.2f %5.0f %5.0f %8.3f %8.2f %8.2f %8.3f %8.3f %8.2f %9.1f\n', ...
        pts{k,1}, pts{k,2}, pts{k,3}, pts{k,4}, rmsv, ymed, yp95, Smean, Sp50, Sp95, Cs);
    fprintf(fid, '%-5s %5.2f %5.0f %5.0f %8.3f %8.2f %8.2f %8.3f %8.3f %8.2f %9.1f\n', ...
        pts{k,1}, pts{k,2}, pts{k,3}, pts{k,4}, rmsv, ymed, yp95, Smean, Sp50, Sp95, Cs);
end
fclose(fid);

fprintf('\n噪声: mean = %g, sigma = %g (pulses/s)\n', 90, 30);
fprintf('S 的饱和值 2*e0 = %g\n', 2*e0);
fprintf('C* = noise_mean / mean(S(y_out))，即"单条权重为 1 的连接恰好贡献一个噪声均值"的 C\n');
fprintf('结果已写入 tests/results/calib_c_scale.txt\n');
end
