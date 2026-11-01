function probe_cycle()
% PROBE_CYCLE  两件事：
%  (1) 复位容差 gain_reset_tol 会不会改动 200 s 的存量结果（应当完全不动）；
%  (2) 打开容差之后，长时程里能不能出现"发作—恢复—再发作"的循环。
%
% 背景：legacy 律的三个 k_*1 在增益触底时被置 0，只有三个增益与静息值
% **精确相等**才复位；浮点下几乎不发生，于是每个节点一生只发一次
% （probe_spontaneous 里 600 s 只有 1 次）。wendling 新增 gain_reset_tol
% （默认 0 = 旧行为），本脚本验证容差 0.05 的效果。
%
% 输出：tests/results/probe_cycle.txt

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); setup_paths();
resdir = fullfile(root, 'tests', 'results');

fid = fopen(fullfile(resdir, 'probe_cycle.txt'), 'w');

%% ---------- (1) 存量结果是否被改动 ----------
fprintf(fid, '# (1) ESF case (3.25,22,10)->(4.25,22,10) @10 s, 200 s\n');
fprintf(fid, '# tol     Ge_end   Gs_end   Gf_end  max|dGe| vs tol=0\n');
c0 = run_case(3.25, 22, 10, 4.25, 200, 0);
for tol = [0.02 0.05 0.1]
    c = run_case(3.25, 22, 10, 4.25, 200, tol);
    dE = max(abs(c.Ge - c0.Ge)); dS = max(abs(c.Gs - c0.Gs)); dF = max(abs(c.Gf - c0.Gf));
    str = sprintf('  %-6.2f %8.3f %8.3f %8.3f  %.2e (E) %.2e (S) %.2e (F)', ...
        tol, c.Ge(end), c.Gs(end), c.Gf(end), dE, dS, dF);
    fprintf('%s\n', str); fprintf(fid, '%s\n', str);
end

%% ---------- (2) 长时程循环 ----------
fprintf(fid, '\n# (2) no ESF, legacy + tol=0.05, T=1200 s\n');
fprintf(fid, '# case            n_burst  first   period   Ge_max  Ge_end  Gs_end  Gf_end\n');
cases = { ...
    'A4.60_B30_G25', 4.60, 30, 25; ...
    'A4.60_B25_G10', 4.60, 25, 10; ...
    'A4.40_B25_G15', 4.40, 25, 15};
for k = 1:size(cases,1)
    c = run_case(cases{k,2}, cases{k,3}, cases{k,4}, [], 1200, 0.05);
    w = round(2*c.fs); nr = floor(numel(c.y)/w);
    rr = zeros(1,nr);
    for i = 1:nr, seg = c.y((i-1)*w+(1:w)); rr(i) = sqrt(mean(seg.^2)); end
    tr = (0:nr-1)*2 + 1;
    base = median(rr(tr <= 20));
    above = rr > max(3*base, 1.0);
    d = diff([0 above 0]);
    nb = sum(d == 1);
    on = tr(d(1:end-1) == 1);
    if numel(on) >= 2, per = median(diff(on)); else, per = NaN; end
    if isempty(on), tf = NaN; else, tf = on(1); end
    str = sprintf('  %-14s %5d %7.1f %8.1f %8.2f %8.2f %8.2f %8.2f', ...
        cases{k,1}, nb, tf, per, max(c.Ge), c.Ge(end), c.Gs(end), c.Gf(end));
    fprintf('%s\n', str); fprintf(fid, '%s\n', str);
end
fclose(fid);
fprintf('\n已写 results/probe_cycle.txt\n');
end

% ----------------------------------------------------------------------
function c = run_case(A, B, G, esf_to, T, tol)
cfg = default_config();
cfg.node_count = 1;
cfg.A = A; cfg.B = B; cfg.G = G;
cfg.model_dynamic = 1;
cfg.gain_law = 'legacy';
cfg.gain_reset_tol = tol;
cfg.simulation_time = T;
cfg.resample_dy = true;
cfg.seed = 20261101;
if isempty(esf_to)
    cfg.esf_enable = false;
else
    cfg.esf_enable = true;
    cfg.esf_time = 10;
    cfg.esf_params = [esf_to, B, G];
end
[y, dy] = simulate_network(cfg);
c.y = y; c.fs = cfg.eeg_fs;
c.Ge = dy(1,:); c.Gs = dy(2,:); c.Gf = dy(3,:);
end
