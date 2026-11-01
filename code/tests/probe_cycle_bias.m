function probe_cycle_bias()
% PROBE_CYCLE_BIAS  用"背景放电基线" gain_bias 让静息点自己爬，检验能否出现
% 发作—崩塌—恢复—再发作的循环。
%
% 依据：静息 (3.25,22,10) 的驱动 x 恒为 0，增益一动不动；而一旦恢复后回到静息，
% 就再也起不来（probe_cycle: 1200 s 只有 1 次）。gain_bias 给 x 加一个常数基线，
% 使 G_e 爬向 rest + k_E1·x0/k_E3 = 3.25 + 200·x0，停在阈值附近 → 反复越界。
%
% 扫 x0 = 0.002/0.003/0.005，另加两个"抗性"参数点作对照（它们阈值更高，应当不动）。
%
% 输出：tests/results/cycle_bias.txt

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); setup_paths();
resdir = fullfile(root, 'tests', 'results');

T = 1200;
cases = { ...
    'SOZ_b2m',  3.25, 22, 10, 0.002; ...
    'SOZ_b3m',  3.25, 22, 10, 0.003; ...
    'SOZ_b5m',  3.25, 22, 10, 0.005; ...
    'PZ_b3m',   3.25, 30, 10, 0.003; ...
    'NIZ_b3m',  3.25, 45, 25, 0.003};

fid = fopen(fullfile(resdir, 'cycle_bias.txt'), 'w');
fprintf(fid, '# no ESF, legacy + reset_tol=0.05, T=%d s\n', T);
fprintf(fid, '# case      n_burst  first   period  dur_med  Ge_max  Ge_min  Ge_end  Gs_end  Gf_end\n');

for k = 1:size(cases,1)
    cfg = default_config();
    cfg.node_count = 1;
    cfg.A = cases{k,2}; cfg.B = cases{k,3}; cfg.G = cases{k,4};
    cfg.gain_bias = cases{k,5};
    cfg.model_dynamic = 1;
    cfg.gain_law = 'legacy';
    cfg.gain_reset_tol = 0.05;
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
    above = rr > max(3*base, 1.0);
    d = diff([0 above 0]);
    nb = sum(d == 1);
    on = tr(d(1:end-1) == 1);
    off = tr(d(2:end) == -1);
    if numel(on) >= 2, per = median(diff(on)); else, per = NaN; end
    if isempty(on), tf = NaN; else, tf = on(1); end
    if numel(on) >= 1 && numel(off) >= 1
        m = min(numel(on), numel(off));
        dur = median(off(1:m) - on(1:m));
    else
        dur = NaN;
    end

    str = sprintf('  %-9s %6d %7.1f %8.1f %8.1f %8.2f %8.2f %8.2f %8.2f %8.2f', ...
        cases{k,1}, nb, tf, per, dur, max(Ge), min(Ge), Ge(end), Gs(end), Gf(end));
    fprintf('%s\n', str); fprintf(fid, '%s\n', str);
end
fclose(fid);
fprintf('\n已写 results/cycle_bias.txt\n');
end
