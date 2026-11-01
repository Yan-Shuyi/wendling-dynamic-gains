function [d, y, Fs_eff] = slope_series(x, Fs, lp)
% SLOPE_SERIES  低通 → 统一重采样到 512 Hz → |dy/dt|（mV/s）。
%
%   [d, y, Fs_eff] = slope_series(x, Fs, lp)
%
%   d       |dy/dt|（mV/s），长度 = numel(y) - 1
%   y       低通并重采样后的信号（与 d 同一时间基准）
%   Fs_eff  重采样后的采样率（512）
%
% 为什么必须在同一个函数里做这两步
%   斜率是带宽依赖的量，逐点差分还强烈依赖采样率：同一个 2 Hz 棘波串在 128 Hz
%   下只有 2.6 个采样点，斜率统计随机跳。参考段与待判信号如果不在同一采样率下
%   算斜率，阈值就会随 Fs 漂移（实测 ref.p99 从 128 Hz 的 190 漂到 1024 Hz 的 70）。
%   两边都走这一个函数，阈值和统计量才落在同一把尺子上。
x = double(x(:));
y = bandlimit(x - mean(x), Fs, lp);
Fsd = 512;
if abs(Fs - Fsd) > 1
    t_src = (0:numel(y)-1)'/Fs;
    t_dst = (0:1/Fsd:t_src(end))';
    y = interp1(t_src, y, t_dst, 'linear', 'extrap');
    Fs = Fsd;
end
Fs_eff = Fs;
d = abs(diff(y)) * Fs;
end
