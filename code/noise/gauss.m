function input_signal = gauss(duration, dt, mu, sigma, method)
% GAUSS  Generate the external afferent drive p(t) as a Gaussian process.
%
%   input_signal = gauss(duration, dt, mu, sigma, method)
%
% 输入
%   duration  时长（s）
%   dt        采样间隔（s），须与仿真步长一致
%   mu        均值，单位 pulses/s；论文用 90
%   sigma     标准差；论文用 30
%   method    分布类型（可选，默认 'normal'）：
%               'normal'     普通高斯，可能出现负值
%               'truncated'  截断到非负
%               'folded'     折叠（取绝对值）
%               'lognormal'  对数正态，恒正
%
% 输出
%   input_signal  1 × (ceil(duration/dt)+1) 的行向量
%
% 做什么
%   生成模型的外部传入驱动。论文正文写的是 p(t) ~ N(90, 30²)，
%   对应 method = 'normal'。
%
% 调用关系
%   上游：gauss_noise（按段调用，本函数只管"一段"怎么抽）
%
% 备注
%   本函数自己不设种子，可复现性由 simulate_network 开头的 rng(cfg.seed) 保证。
%   'lognormal' 依赖 lognrnd（Statistics and Machine Learning Toolbox）；
%   论文三张图都用 'normal'，不需要该工具箱。
%   注意 'normal' 允许负值：mu = 90、sigma = 30 时约 0.13% 的样本为负，
%   即"负的传入脉冲密度"。论文用 'normal' 是照原文做法，但严格说 p(t) 应恒正，
%   若要较真就换 'truncated'（会改变图，需重跑）。
    if nargin < 5
        method = 'normal'; % 推荐
    end

    n_samples = ceil(duration / dt) + 1;

    switch method
        case 'normal'
            input_signal = mu + sigma * randn(1, n_samples);

        case 'truncated'
            input_signal = mu + sigma * randn(1, n_samples);
            input_signal(input_signal < 0) = 0;

        case 'folded'
            input_signal = mu + abs(sigma * randn(1, n_samples));

        case 'lognormal'
            sigma2 = log(1 + (sigma/mu)^2);
            mu2 = log(mu) - sigma2/2;
            input_signal = lognrnd(mu2, sqrt(sigma2), 1, n_samples);
    end
end