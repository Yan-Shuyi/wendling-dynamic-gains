function y = bandlimit(x, Fs, lp)
% BANDLIMIT  零相位低通（4 阶 Butterworth，filtfilt）；lp<=0 或 lp>Nyquist 时原样返回。
%
%   y = bandlimit(x, Fs, lp)
%
% 为什么需要它
%   斜率（mV/s）是带宽依赖的量：同一个信号，采样率越高逐点差分越大。
%   算出任何"斜率类"判据之前先把高频裁到生理带宽内，判据才和 Fs 无关，
%   也才不会被参考段里的高频噪声绑架（见 detect_spikes 的说明）。
x = double(x(:));
if lp <= 0 || lp >= Fs/2
    y = x;
    return;
end
[b, a] = butter(4, lp/(Fs/2));
y = filtfilt(b, a, x);
end
