function test_classifier_archetypes()
% TEST_CLASSIFIER_ARCHETYPES  用"应判成什么"已知的合成原型，对照 v1（旧树）与 v2（文献判据）。
%
%   运行：run_checks   （或在任意目录执行 addpath('<code>'); setup_paths; test_classifier_archetypes）
%
% 为什么要有这个测试
%   旧树从没做过任何验证。这里每条原型的期望标签由波形本身决定（不是手标），
%   任一不符就是判据有结构问题。重点看三条：
%     - 1 Hz 大幅正弦（应 ICTAL 慢节律）在 v1 里会被判成"正常背景"
%     - 40 Hz 无偏置快活动（应 ONSET）在 v1 里会被判成"持续棘波"
%     - 同一信号换采样率，v1 的类别会跳；v2 不应跳
%
% 输出：results/classifier_archetypes.txt
here = fileparts(mfilename('fullpath'));
outdir = fullfile(here, 'results');
if ~exist(outdir, 'dir'), mkdir(outdir); end
setup_paths();
fid = fopen(fullfile(outdir, 'classifier_archetypes.txt'), 'w');
p1 = classify_params();          % 旧
p2 = classify_params_v2();       % 新
Fs = 256;

%% 参考段：低幅背景噪声（模拟发作间期基线）
rng(7);
base = 0.8*randn(round(p2.win_sec*Fs),1) + 0.5*sin(2*pi*1*(0:1/Fs:p2.win_sec-1/Fs)');
ref = make_reference(base, Fs, p2);

%% 原型列表：{名字, 信号, 期望标签(v2), 期望标签(v1)}
T = 10; t = (0:1/Fs:T-1/Fs)';
c = {};
c(end+1,:) = {'1 Hz 正弦 amp=6, dc=0',        6*sin(2*pi*1*t),        'uncertain-offband', 'Type1'};
c(end+1,:) = {'4 Hz 正弦 amp=8, dc=12',       12+8*sin(2*pi*4*t),     'ICTAL-theta',    'Type4'};
c(end+1,:) = {'9 Hz 正弦 amp=5, dc=10',       10+5*sin(2*pi*9*t),     'ICTAL-alpha',    'Type4'};
c(end+1,:) = {'25 Hz 快活动 amp=2, dc=0',     2*sin(2*pi*25*t),       'ONSET-lvf',      'Type1'};
c(end+1,:) = {'35 Hz 快活动 amp=2, dc=8',     8+2*sin(2*pi*35*t),     'ONSET-lvf',      'Type5'};
c(end+1,:) = {'2 Hz 棘波串 amp=8, dc=0',      spike_train(t,2,0.02,8),'preONSET-repetitive','Type3/2'};
c(end+1,:) = {'0.3 Hz 棘波串 amp=8, dc=0',    spike_train(t,0.3,0.02,8),'preONSET-sporadic','Type2/3'};
c(end+1,:) = {'白噪声 sigma=0.8, dc=0',       0.8*randn(size(t)),     'intICTAL',       'Type2/3'};
c(end+1,:) = {'常数 10 mV',                   10*ones(size(t)),       'intICTAL',       'Type6'};

fprintf(fid, '原型对照（win=10 s, Fs=%d，参考段=低幅背景）\n', Fs);
fprintf(fid, '%-30s | %-22s | %-10s | %-10s\n', '原型', '期望(v2)', 'v2 判定', 'v1 判定');
fprintf(fid, '%s\n', repmat('-',1,82));
for k = 1:size(c,1)
    nm = c{k,1}; x = c{k,2}; exp2 = c{k,3}; exp1 = c{k,4};
    [lab, sub] = classify_v2(x, Fs, p2, ref);
    if isempty(sub), s2 = lab; else, s2 = [lab '-' sub]; end
    % 旧树：5 s 窗 + 它自己的阈值
    m = min(numel(x), round(5*Fs));
    f1 = extract_features(x(1:m), p1);
    t1 = classify_single(f1, p1);
    fprintf(fid, '%-30s | %-22s | %-10s | Type %d\n', nm, exp2, s2, t1);
end

%% 采样率不变性
fprintf(fid, '\n采样率不变性（同一 2 Hz 棘波串，只改 Fs）\n');
for Fs2 = [128 256 512 1024]
    tt = (0:1/Fs2:10-1/Fs2)';
    x = spike_train(tt, 2, 0.02, 8);            % 同样的 2 Hz 棘波串（20 个）
    rng(3);                                     % 参考段噪声每次同分布
    ref2 = make_reference(0.8*randn(round(p2.win_sec*Fs2),1), Fs2, p2);
    [lab, sub, f2] = classify_v2(x, Fs2, p2, ref2);
    f1 = extract_features(x, setfield(p1,'Fs',Fs2)); %#ok<SFLD>
    fprintf(fid, '  Fs=%5d | v2: %-11s %-11s spk=%.2fHz(真值2.0) refp99=%6.1f | v1: Type %d (diff_max=%.2f)\n', ...
        Fs2, lab, sub, f2.spk_rate, ref2.slope_p99, classify_single(f1, p1), f1.diff_max);
end
fclose(fid);
fprintf('已写出：%s\n', fullfile(outdir, 'classifier_archetypes.txt'));
end

function x = spike_train(t, rate, width, amp)
x = zeros(size(t));
for c0 = 0.1:(1/rate):(t(end)-0.1)
    x = x + amp*exp(-((t-c0)/width).^2);
end
end
