function test_onset_alignment()
% TEST_ONSET_ALIGNMENT  onset 定义敏感性 + 中位数/均值两种汇总的差异（原 .tmp_d13/d14c）。
%
%   test_onset_alignment
%
% 回答两个问题
%   1. 把 onset 判据从 Gs<10 换成中点 12.75 或 8/15，结论（四期对齐后中位数一致）会不会变？
%   2. 同一个量用中位数+IQR 报 vs 用 mean±SD 报，差多少？为什么选前者？
%
% 输出：results/onset_alignment.txt
here = fileparts(mfilename('fullpath'));
outdir = fullfile(here, 'results');
if ~exist(outdir, 'dir'), mkdir(outdir); end
setup_paths();

fid = fopen(fullfile(outdir, 'onset_alignment.txt'), 'w');
offs = [-10 0 10 20]; nm = {'-10','0','+10','+20'};
thrs = [8 10 12.75 15];
fprintf(fid, 'onset 定义敏感性：Gs 首次跌破 thr（静息 22，下限约 3.5，中点 12.75）\n');

for ge = [4.25 4.00]
    cfg = default_config();
    cfg.node_count = 20; cfg.esf_enable = true; cfg.esf_time = 10;
    cfg.esf_params = [ge 22 10]; cfg.esf_jitter = 0.05;
    [~, dy, cfg] = simulate_network(cfg);
    fs = cfg.dy_fs; n = size(dy, 2);
    fprintf(fid, '\n=== esf Ge = %.2f ===\n', ge);

    for thr = thrs
        T = NaN(cfg.node_count, 1);
        for i = 1:cfg.node_count
            k = find(dy((i-1)*3+2, :) < thr, 1, 'first');
            if ~isempty(k), T(i) = (k-1)/fs; end
        end
        ok = ~isnan(T);
        fprintf(fid, ' thr=%5.2f  检出 %2d/%d  onset 中位 %.1f s (%.0f-%.0f)\n', ...
            thr, sum(ok), cfg.node_count, median(T(ok)), min(T(ok)), max(T(ok)));
        if thr ~= 12.75, continue; end

        % 中点阈值下给出完整对照表
        for g = 1:3
            gn = {'Ge','Gs','Gf'};
            fprintf(fid, '   %s  off     n    median   IQR     mean    SD\n', gn{g});
            for q = 1:numel(offs)
                v = [];
                for i = find(ok)'
                    k = round((T(i)+offs(q))*fs)+1;
                    if k >= 1 && k <= n, v(end+1) = dy((i-1)*3+g, k); end %#ok<AGROW>
                end
                fprintf(fid, '        %3s  %4d  %7.2f  %6.2f  %7.2f  %5.2f\n', nm{q}, numel(v), ...
                    median(v), prctile(v,75)-prctile(v,25), mean(v), std(v));
            end
        end
    end
end
fclose(fid);
fprintf('已写出：%s\n', fullfile(outdir, 'onset_alignment.txt'));
end
