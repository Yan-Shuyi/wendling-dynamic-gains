function outputs_fs = resample_fs(outputs, originalFs, targetFs, use_parallel)
% RESAMPLE_FS  Resample a multi-channel matrix to a target sampling rate.
%
%   outputs_fs = resample_fs(outputs, originalFs, targetFs, use_parallel)
%
% 输入
%   outputs       原始数据，[通道数 × 采样点数]
%   originalFs    原始采样率（Hz），仿真里是 1/dt = 256
%   targetFs      目标采样率（Hz）；EEG 用 256（不降采样），增益轨迹用 dy_fs
%   use_parallel  是否用 parfor 并行（可选，默认 false）
%
% 输出
%   outputs_fs    重采样后的 [通道数 × 新采样点数]
%
% 做什么
%   逐通道调用 MATLAB 内置 resample（带抗混叠 FIR 滤波），避免直接抽点造成
%   频谱混叠。通道维度固定为第 1 维。
%
% 调用关系
%   上游：simulate_network（对 EEG 与增益轨迹各调一次）
%
% 备注
%   100 节点 × 200 s 的 fig3 里，本函数是耗时的主要来源之一。
    if nargin < 4
        use_parallel = false;
    end
    
    n_channels = size(outputs, 1);  % 修正：第一维是通道数
    
    if ~use_parallel
        % 串行处理：直接沿第二维（时间维）重采样
        outputs_fs = resample(outputs, targetFs, originalFs, 'Dimension', 2);
    else
        % 并行处理
        if ~license('test', 'Distrib_Computing_Toolbox')
            warning('Parallel Computing Toolbox 不可用，将使用串行处理。');
            outputs_fs = resample(outputs, targetFs, originalFs, 'Dimension', 2);
            return;
        end
        
        try
            if isempty(gcp('nocreate'))
                parpool('local');
            end
        catch ME
            warning('无法启动并行池，将使用串行处理');
            outputs_fs = resample(outputs, targetFs, originalFs, 'Dimension', 2);
            return;
        end
        
        % 预分配输出数组
        temp_resampled = resample(outputs(1,:), targetFs, originalFs);
        n_samples_new = length(temp_resampled);
        outputs_fs = zeros(n_channels, n_samples_new);
        
        % 并行处理每个通道
        parfor i = 1:n_channels
            outputs_fs(i, :) = resample(outputs(i,:), targetFs, originalFs);
        end
    end
end