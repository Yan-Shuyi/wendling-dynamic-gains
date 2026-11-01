function ref = make_reference(x, Fs, p)
% MAKE_REFERENCE  从参考段（同一仿真的发作间期 / 或纯噪声实现）算归一化基准。
%
%   ref = make_reference(x, Fs, p)
%
% 为什么需要它
%   v1 的阈值是绝对 mV，换噪声均值、换增益律、换状态机 k 都要重调，且不报错。
%   这里把所有幅度类判据改成"相对参考段"的比值，一次标定、到处适用。
%
% 输出
%   ref.rms        参考段 RMS（mV）
%   ref.range      参考段 p95-p5（mV），比 max-min 抗野点
%   ref.slope_p99  参考段 |dx/dt| 的 99 百分位（mV/s），棘波阈值由它定
%   ref.psd        参考段 1–90 Hz 平均功率谱（用于将来做谱距离）
%   ref.f          对应频率轴（Hz）
x = x(:);
if p.hp_cut > 0
    x = highpass(x, p.hp_cut, Fs);
else
    x = x - mean(x);
end

ref.rms = sqrt(mean(x.^2));
ref.range = prctile(x, 95) - prctile(x, 5);

% 斜率先裁到生理带宽（spk_lp）并统一到 512 Hz，与待判信号走同一个函数
d = slope_series(x, Fs, p.spk_lp);   % mV/s
ref.slope_p99 = prctile(d, 99);
if ref.slope_p99 <= 0
    ref.slope_p99 = eps;
end

win = min(numel(x), round(2*Fs));
nfft = max(512, 2^nextpow2(win));
[pxx, f] = pwelch(x, hanning(win), floor(win/2), nfft, Fs);
keep = f >= 1 & f <= 90;
ref.psd = pxx(keep);
ref.f   = f(keep);
end
