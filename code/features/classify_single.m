% Start [输入5s信号]
%   │
%   ├─── 1. 计算 Mean (均值)
%   │
%   └─── Mean >= 5mV ? ──YES──▶ [高幅组 4,5,6]
%                                 │
%                                 ├─► 高频占比高? ──YES──▶ Type 5
%                                 │
%                                 ├─► 带宽窄? ───────YES──▶ Type 6
%                                 │
%                                 └─► 否则 ───────────────▶ Type 4
%                               (绝不检查Min/负值)
%   │
%  NO
%   │
%   ▼
% [低幅组 1,2,3] (这里才允许负值存在)
%   │
%   ├─── 2. 计算 Range (极差)
%   │
%   └─── Range < 5mV ? ──YES──▶ Type 1 (平稳背景)
%         │
%        NO
%         │
%         ▼
%   ├─── 3. 计算 Diff_Max (差分)
%         │
%         └─── Diff_Max > Thresh ? ──NO──▶ Type 1 (虽有波动但无陡变)
%               │
%              YES
%               │
%               ▼
%         ├─── 4. 计算 Spectrum Peaks (频谱尖峰)
%               │
%               ├─── Peaks <= 2 ──▶ Type 3 (持续棘波)
%               │
%               └─── Peaks > 2 ───▶ Type 2 (零星棘波)
% ==============================
% plaintext: 单窗分类逻辑（严格执行最终规则）
% 逻辑链：
% 1. 先看均值：>=5mV→高幅组[4,5,6]（绝不检查最小值/负值）
% 2. 均值<5mV→低幅组[1,2,3]（此处才允许讨论负值）
%   2.1 极差<2mV→Type1（平稳背景）
%   2.2 极差>=2mV→查差分：diff_max>3*std(diff)→有棘波
%     2.2.1 频谱尖峰<=2→Type3（持续棘波，纯净）
%     2.2.2 频谱尖峰>2→Type2（零星棘波，杂乱）
%     2.2.3 无陡变→Type1（缓慢漂移）
% ==============================
function type_label = classify_single(features, params)
% CLASSIFY_SINGLE  Judge one EEG window and return its state type (1-6).
%
%   [type_id, feat, branch] = classify_single(x, Fs, params)
%
% 输入
%   x       一个滑窗内的 EEG 片段（1 × n，mV）
%   Fs      采样率（Hz）
%   params  classify_params() 的阈值集合
%
% 输出
%   type_id  1..6 的类别编号；判不出来时返回 NaN
%   feat     extract_features 算出的全部特征（便于事后核对为什么判成这一类）
%   branch   走了哪条枝（'high' / 'low'），调试用
%
% 做什么（就是论文 Fig. 2 那棵树）
%   第一层：均值 >= thresh_mean(5 mV) → 高幅组 {4,5,6}，否则低幅组 {1,2,3}
%   高幅枝：高频能量占比高 → 5；90% 带宽窄 → 6；否则 → 4
%   低幅枝：极差小 → 1；最大一阶差分小 / 显著峰少 → 2；否则 → 3
%   所有阈值都来自 classify_params.m，论文 Fig. 2 中的数字与之一一对应。
%
% 调用关系
%   上游：classify_eeg_timeseries（逐窗调用）
%   下游：extract_features
%
% 备注
%   文件顶部原有的 ASCII 流程图与代码已不一致（图上写 5 mV 两级、实际多出
%   带宽与峰数判据），以本头注释和代码为准。
type_label = 1; % 默认Type1

%% ----------------------
% 1. 一级分类：均值一刀切（核心规则，无交叉）
% ----------------------
if features.mean_val >= params.thresh_mean_high
    % ----------------------
    % 高幅组[4,5,6]：绝不检查最小值/负值，仅用频域特征
    % ----------------------
    if features.high_freq_ratio > params.thresh_high_freq_ratio
        type_label = 5; % 低幅快活动（高频主导）
    else
        if features.spec_bw < params.thresh_spec_bw
            type_label = 6; % 高幅准正弦（带宽极窄）
        else
            type_label = 4; % 慢节律（带宽较宽）
        end
    end
    return; % 高幅组逻辑终结，绝不进入低幅组判断
end

%% ----------------------
% 2. 低幅组[1,2,3]：此处才允许讨论负值、极差、差分
% ----------------------
% 2.1 极差极小→平稳背景（优先级最高）
if features.range_val < params.thresh_range_123
    type_label = 1;
    return;
end

% 2.2 极差大，检查是否有棘波陡变
if features.diff_max <= params.diff_max
    % 无陡变：即使是大幅波动，也只是缓慢漂移，归为背景
    type_label = 1;
    return;
end

% 2.3 有陡变（棘波），按频谱纯净度细分2和3
if features.peak_num <= params.thresh_peak_num
    type_label = 3; % 持续棘波：频谱纯净，峰少
else
    type_label = 2; % 零星棘波：频谱杂乱，峰多
end
end