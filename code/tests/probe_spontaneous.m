function probe_spontaneous()
% PROBE_SPONTANEOUS  在"驱动占空比很小"的点上，慢动力学能不能自己把节点推入发作，
% 以及发作之后会不会有第二次。
%
% 依据 probe_rest_x 的占空比表：静息 (3.25,22,10) 的 duty = 0（x 恒为 0，增益不动），
% 而 A=4.00/B=25/G=25 的 duty = 0.052（约 5% 的时间 x>0），A=4.60/B=30/G=25 是 0.28。
% 本脚本在这三个点上跑 600 s、不给 ESF，记录：
%   首次放电时刻、放电次数、Ge 的极值与末值、以及相邻两次放电的间隔。
%
% 输出：tests/results/spontaneous.txt

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); setup_paths();
resdir = fullfile(root, 'tests', 'results');

T = 600;
cases = { ...
    'A4.00_B25_G25', 4.00, 25, 25; ...
    'A4.60_B30_G25', 4.60, 30, 25; ...
    'A4.60_B25_G10', 4.60, 25, 10};

fid = fopen(fullfile(resdir, 'spontaneous.txt'), 'w');
fprintf(fid, '# no ESF, legacy gain law, T=%d s\n', T);
fprintf(fid, '# case            Ge_min Ge_max Ge_end Gs_max Gs_min Gf_max  n_burst first  period\n');

for k = 1:size(cases, 1)
    cfg = default_config();
    cfg.node_count = 1;
    cfg.A = cases{k,2}; cfg.B = cases{k,3}; cfg.G = cases{k,4};
    cfg.model_dynamic = 1;
    cfg.gain_law = 'legacy';
    cfg.esf_enable = false;
    cfg.simulation_time = T;
    cfg.resample_dy = true;
    cfg.seed = 20261101;

    [y, dy] = simulate_network(cfg);
    fs = cfg.eeg_fs;
    Ge = dy(1,:); Gs = dy(2,:); Gf = dy(3,:);

    w = round(2*fs); nr = floor(numel(y)/w);
    rr = zeros(1,nr);
    for i = 1:nr, seg = y((i-1)*w+(1:w)); rr(i) = sqrt(mean(seg.^2)); end
    tr = (0:nr-1)*2 + 1;
    base = median(rr(tr <= 20));
    above = rr > max(3*base, 1.0);       % 同时要求绝对幅度 > 1 mV，避开噪声地板
    d = diff([0 above 0]);
    nb = sum(d == 1);
    on = tr(d(1:end-1) == 1);
    if numel(on) >= 2, per = median(diff(on)); else, per = NaN; end
    if isempty(on), tf = NaN; else, tf = on(1); end

    str = sprintf('%-14s %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f  %7d %6.1f %6.1f', ...
        cases{k,1}, min(Ge), max(Ge), Ge(end), max(Gs), min(Gs), max(Gf), nb, tf, per);
    fprintf('%s\n', str);
    fprintf(fid, '%s\n', str);
end
fclose(fid);
fprintf('\n已写 results/spontaneous.txt\n');
end
