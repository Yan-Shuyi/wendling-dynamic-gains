function run_checks()
% RUN_CHECKS  一次跑完所有"不画图、只出数"的验证脚本，结果写到 code/tests/results/。
%
%   用法（任意工作目录下都能跑）
%     addpath('C:/Users/26442/WorkBuddy/Claw/code'); run_checks
%     或先 cd 到 code/ 再 run_checks
%
% 跑什么
%   1. test_classifier_archetypes  v1/v2 分期判据的原型断言 + 采样率不变性
%   2. test_classifier_boundaries  判据边界过渡 / 噪声鲁棒性 / 独立描述子交叉检验
%                                   + 适用域（多窄的瞬变会被误当快活动）
%   3. test_onset_alignment        onset 定义敏感性 + 中位数/均值汇总对照
%
% 为什么单独做成一个入口
%   之前的脚本散在临时目录里，用户粘贴路径容易出错（反斜杠被吞、文件夹名带点）。
%   放在 code/tests 下、由 setup_paths 挂进路径，就只需要一条 addpath。
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));      % code/
setup_paths();

fprintf('\n===== 1/3  分期判据原型检验 =====\n');
try
    test_classifier_archetypes();
catch err
    fprintf(2, '失败：%s\n', err.message);
end

fprintf('\n===== 2/3  判据边界与鲁棒性 =====\n');
try
    test_classifier_boundaries();
catch err
    fprintf(2, '失败：%s\n', err.message);
end

fprintf('\n===== 3/3  onset 对齐与汇总统计 =====\n');
try
    test_onset_alignment();
catch err
    fprintf(2, '失败：%s\n', err.message);
end

fprintf('\n结果目录：%s\n', fullfile(here, 'results'));
end
