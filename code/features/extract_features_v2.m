function f = extract_features_v2(x, Fs, p)
% EXTRACT_FEATURES_V2  分期判据所需的特征：频带功率、主频、棘波率、去极化。
%
%   f = extract_features_v2(x, Fs, p)
%
% 与 v1 的差别
%   - 先减均值 / 高通，再做一切（v1 用均值当幅度，对称大幅振荡会掉进低幅组）
%   - 棘波用斜率（mV/s）+ 不应期计数，得到**棘波率**，不再用"谱峰个数"（v1 测的是谐波）
%   - 幅度用 rms / p95-p5，不用 max-min
%   - 直流偏置单独输出 f.dc，只作描述、不参与分期
x = double(x(:));
if p.hp_cut > 0
    xf = highpass(x, p.hp_cut, Fs);
else
    xf = x - mean(x);
end
f.dc  = mean(x);              % 去极化程度（只描述，不判类）
f.n   = numel(xf);

%% 时域
f.rms   = sqrt(mean(xf.^2));
f.range = prctile(xf, 95) - prctile(xf, 5);
% 斜率与棘波一样先裁到 spk_lp，保证与 Fs 无关、也与参考段高频无关
[d, xb] = slope_series(xf, Fs, p.spk_lp);
f.crest = max(abs(xb)) / max(sqrt(mean(xb.^2)), eps);   % 波峰因数：节律 ≈1.4，棘波 ≥3
f.slope_max = max(d);
f.slope_p99 = prctile(d, 99);

%% 频域
win = min(numel(xf), round(2*Fs));
nfft = max(512, 2^nextpow2(win));
[pxx, fr] = pwelch(xf, hanning(win), floor(win/2), nfft, Fs);
keep = fr >= 1 & fr <= 90;
pxx = pxx(keep); fr = fr(keep);
tot = sum(pxx);
if tot <= 0, tot = eps; end
f.f = fr; f.psd = pxx;

% 分带功率（绝对 + 相对）
nb = size(p.bands, 1);
f.abs = zeros(nb,1); f.rel = zeros(nb,1);
for b = 1:nb
    m = fr >= p.bands{b,2} & fr < p.bands{b,3};
    f.abs(b) = sum(pxx(m));
    f.rel(b) = f.abs(b) / tot;
    f.(p.bands{b,1}) = f.rel(b);       % 便于按名字取：f.theta 等
end
% 两个判据带（可能与 bands 不重合，单独算）
mF = fr >= p.fast_band(1) & fr <= p.fast_band(2);
mS = fr >= p.slow_band(1) & fr <= p.slow_band(2);
f.rel_fast = sum(pxx(mF)) / tot;       % ONSET 判据
f.rel_slow = sum(pxx(mS)) / tot;       % ICTAL 判据

[~, ip] = max(pxx);
f.f_peak = fr(ip);

% 90% 能量带宽，用 f10/f90 两端，不用"减 1 Hz"那种写法
c = cumsum(pxx)/tot;
f.f10 = fr(find(c >= 0.10, 1, 'first'));
f.f90 = fr(find(c >= 0.90, 1, 'first'));
if isempty(f.f10), f.f10 = fr(1); end
if isempty(f.f90), f.f90 = fr(end); end
f.bw90 = f.f90 - f.f10;

end
