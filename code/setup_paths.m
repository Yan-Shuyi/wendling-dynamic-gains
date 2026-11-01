function [root, figdir] = setup_paths()
%SETUP_PATHS  Add the active code folders to the MATLAB path.
%
%   [root, figdir] = setup_paths()
%
% 输入
%   无
%
% 输出
%   root     项目 code/ 目录的绝对路径（本文件所在目录）
%   figdir   论文插图目录（项目根下的 figures/），所有 run_*.m 都往这里写图
%
% 做什么
%   1. 把 core / model / network / noise / features / plotting / util 加入搜索路径；
%   2. 把 code/ 本身也加入，方便直接调用 run_* 脚本；
%   3. legacy/ 与 archive/ 故意不加路径——那是已废弃分支，加了容易误调到同名函数。
%
% 说明
%   每次运行 run_fig*.m 都会先调本函数，所以在任意工作目录下执行都能找到函数，
%   不要求你先手动 cd 到 code/。

root = fileparts(mfilename('fullpath'));

dirs = {'core', 'model', 'network', 'noise', 'features', 'plotting', 'util', 'tests'};
% tests/ 只放验证脚本（如 test_classifier_archetypes），不参与出图；缺失时跳过
for k = 1:numel(dirs)
    d = fullfile(root, dirs{k});
    if exist(d, 'dir')
        addpath(d);
    end
end
addpath(root);

% 论文插图统一目录：<项目根>/figures
figdir = fullfile(root, '..', 'figures');
if ~exist(figdir, 'dir')
    mkdir(figdir);
end

end
