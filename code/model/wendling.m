classdef wendling < base
% WENDLING  Wendling 四群体神经质量模型（含可选的动态增益 ESF）
%
%   m = wendling(G_e, G_s, G_f, dt, model_dynamic)
%
% 输入
%   G_e, G_s, G_f   三个突触增益（mV）：兴奋 / 慢抑制 / 快抑制
%   dt              积分步长（s）
%   model_dynamic   0 = 增益当常数；1 = 增益由 ESF 迟滞律驱动成状态变量
%
% 输出（作为对象属性使用）
%   output          当前时刻的模型输出 y_out = y_p - y_s - y_f（mV），即模拟 EEG
%   state_size = 10 十个状态量 = 五组 (y, x)，对应原文 y0..y9
%   dy_ESF          本次调用的 dG_e/dt, dG_s/dt, dG_f/dt（静态模式下恒为 0）
%
% 做什么
%   方程与 Wendling et al. (2002, Eur. J. Neurosci. 15:1499) 的式(1)逐项对应。
%   原文是 **四个神经元子群**（PC / EIN / SIN / FIN），但状态变量有 10 个，
%   因为慢抑制通路写了两份结构完全相同的二阶方程。状态编号（对应原文 y0..y9）：
%
%       idx 1/6    y0/y5   PC    锥体细胞，模型输出从这里取
%       idx 2/7    y1/y6   EIN   兴奋性中间神经元，外部驱动 p(t) 只进这里
%       idx 3/8    y2/y7   SIN   慢树突抑制，经 C4 反馈到 PC
%       idx 4/9    y3/y8   FIN   快胞体抑制，经 C7 反馈到 PC
%       idx 5/10   y4/y9   SIN'  慢抑制的第二个拷贝
%
%   关于 y4（最容易记错的一处）
%   y4 不是第五个群体。它与 y2 完全同构：同样由 PC 以 C3 驱动、增益同样是
%   (G_s, w_s)、动力学形式一模一样；唯一区别是 y2 的输出进 PC 的 S 函数，
%   而 y4 的输出进 FIN 的 S 函数内部（作为减项，权重 C6）。原文用 Fig.1b 里
%   "GABA_A,slow → GABA_A,fast" 那条箭头表示它，可以理解为"慢抑制中间神经元
%   投向快抑制中间神经元的那一份活动"，是为了让快抑制被慢抑制去抑制（disin-
%   hibition）而引入的中继状态。
%   因此论文正文必须写 four populations / ten state variables；写成
%   five populations 或 eight state variables 都与代码不一致。
%
%   每一组二阶方程 = PSP 核 + 静态 S 型非线性：
%       z = S(v) = 2 e0 / (1 + exp(r (v0 - v)))
%   导数在 derivatives() 里算，积分用四阶 Runge-Kutta（dynamics()）。
%
% 关于耦合常数的命名（重要，容易看晕）
%   代码里的 C_xx 是照搬 Wendling 原文 C1..C7 的，字母顺序并不严格表示
%   "从谁到谁"。实际作用看 derivatives() 里怎么用：
%     C_ep=135    y_p 进入 EIN 的 S 函数之前的突触数        （原文 C1）
%     C_pe=108    上式算出的发放率 z2 反馈进 EIN 方程的权重  （原文 C2 = 0.8·C1）
%     C_sp=33.75  y_p 进入 SIN 的 S 函数之前的突触数        （原文 C3）
%     C_ps=33.75  SIN 发放率反馈权重                        （原文 C4）
%     C_fp=40.5   y_p 进入 FIN 的 S 函数之前的突触数        （原文 C5）
%     C_pf=13.5   y4 在 FIN 的 S 函数内部的抑制权重          （原文 C6）
%     C_fs=108    FIN 发放率反馈权重                        （原文 C7）
%   注意每个中间神经元方程里同时出现两个常数（进 S 前一个、乘发放率一个），
%   PC 方程里一个常数都没有，这是原文本身的结构，不是笔误：
%       y5' = A a S(y1 - y2 - y3)          ← PC，S 内无 C
%       y6' = A a {p(t) + C2 S(C1 y0)}     ← EIN，S 内 C1、S 外 C2
%       y7' = B b C4 S(C3 y0)              ← SIN，S 内 C3、S 外 C4
%       y8' = G g C7 S(C5 y0 - C6 y4)      ← FIN，S 内 C5/C6、S 外 C7
%       y9' = B b S(C3 y0)                 ← SIN'，S 内 C3、S 外无 C
%   论文正文若只写 Σ C_ij S(·)（单个常数），与代码和原文都不一致。
%
% 调用关系
%   上游：simulate_network（逐节点构造）
%   下游：base（父类）
%
% 属性与参数
    properties
        % 时间步长
        dt = 1/256;  % % 采样间隔 (s)，第二篇文献使用256Hz

        % 增益复位的容差（mV）。0 = 旧行为：要求三个增益与静息值**精确相等**
        % 才恢复 k_*1，而浮点下几乎永不成立，于是 legacy 律下每个节点一生只发一次。
        % 取 0.05 时，增益回到静息值附近 0.05 mV 内即复位，长时程仿真里可以出现
        % 第二次、第三次发作（"恢复后等待下一次"）。默认 0，保证 200 s 的存量图不变。
        gain_reset_tol = 0;

        % 背景放电基线 (mV)：加在增益驱动 x 上，见 derivatives 里的注释。
        % 0 = 旧行为（存量图都用 0）；>0 让静息点也能缓慢爬升，用于长时程循环仿真。
        gain_bias = 0;

        % 增益参数
        G_e = 3.25;      % 兴奋性突触增益 (mV)
        G_s = 22;        % 慢抑制性突触增益 (mV)
        G_f = 10;        % 快抑制性突触增益 (mV)

        % 时间常数倒数
        w_e = 100;       % 兴奋性时间常数倒数 (1/s)
        w_s = 50;        % 慢抑制性时间常数倒数 (1/s)
        w_f = 500;       % 快抑制性时间常数倒数 (1/s)

        % 连接常数
        % 注意：字母顺序不严格表示"从谁到谁"，实际作用看 derivatives()。
        % 每条通路有两个常数：一个在 S 函数之前乘 y_p，一个乘算出来的发放率。
        C_ep = 135;          % y_p → EIN 的 S 函数前的突触数          (原文 C1)
        C_pe = 135*0.8;      % 上式发放率反馈进 EIN 方程的权重        (原文 C2)
        C_sp = 135*0.25;     % y_p → SIN 的 S 函数前的突触数          (原文 C3)
        C_ps = 135*0.25;     % SIN 发放率反馈权重                      (原文 C4)
        C_fp = 135*0.3;      % y_p → FIN 的 S 函数前的突触数          (原文 C5)
        C_pf = 135*0.1;      % y4 在 FIN 的 S 函数内部的抑制权重        (原文 C6)
        C_fs = 135*0.8;      % FIN 发放率反馈权重                      (原文 C7)

        % Sigmoid函数参数
        e0 = 2.5;      % 最大放电率 (脉冲/秒)
        v0 = 6;        % 阈值电位 (mV)
        r  = 0.56;     % 斜率参数 (1/mV)

        % ===== 动态参数状态索引（新增3个状态）=====
        if_dy = 0;     % 是否使用动态变量
        dy_size = 0;   % 动态变量维度
        dy_ESF;        % 动态变量
        y_thr = 0.03;  % 变动阈值
        location = [ 'E'; 'S'; 'F'];
        dy_rest  = [3.25;  22;  10]; % 静息态
        dy_min   = [   3;   3;   5]; % 变量最小值
        dy_max   = [   9;  31;  18]; % 变量最大值
        dy_thr   = [ 8.5;  30;  16]; % 上升阈值
        dy_inv   = [   0;   0;   0]; % 反转标志
        x_in     = [   2;  20;   5]; % x 驱动上升系数
        x_de     = [   2;  20;   2]; % x 驱动失活系数
        dy_re    = [0.01;0.01;0.01]; % 恢复系数

        % ===== 增益动力学形式（D13 路径 B）=====
        % 'legacy'    旧的分段线性 + 三个独立状态机（默认，保证三张存量图不变）
        % 'bistable'  统一的三次型双稳方程，三个增益同一个式子，见 gain_bistable.m
        gain_law = 'legacy';

        % 双稳方程参数（gain_law = 'bistable' 时生效）。每行一个增益，顺序 E/S/F。
        % Glo/span 决定 φ 与增益的换算：G = Glo + span·φ
        %   E：φ=0 对应 3 mV，span +6  → φ 大 = 兴奋性增益高（病理方向向上）
        %   S：φ=0 对应 31 mV，span −28 → φ 大 = 慢抑制增益低（病理方向向下）
        %   F：φ=0 对应 18 mV，span −13 → φ 大 = 快抑制增益低（病理方向向下）
        % phi0 取各增益静息值对应的 φ：(3.25, 22, 10) → (0.042, 0.321, 0.615)
        bs = struct( ...
            'Glo',     [3;  31;  18], ...
            'span',    [6; -28; -13], ...
            'phi0',    [0.042; 0.321; 0.615], ...
            'phi_thr', [0.49; 0.77; 0.92], ...
            'A',       [6;  6;  6], ...
            'B',       [0.045; 0.045; 0.030], ...
            'hc',      [0.10; 0.10; 0.12], ...
            'C0',      [0.005; 0.020; 0.015], ...
            'C1',      [0;  0;  0]);

        % ===== E的动力学参数（按你的测试结果调）=====
        E_rest;    % 静息态
        E_min;     % 彻底失活
        E_max;     % 饱和上限
        E_thr;     % E上升阈值
        E_inv;     % 反转
        k_E1;      % x驱动上升系数
        k_E2;      % x驱动失活系数
        k_E3;      % 恢复系数

        % ===== S的动力学参数 =====
        S_rest;    % 静息态
        S_min;     % 彻底失活
        S_max;     % 代偿最大幅值
        S_thr;     % S上升阈值
        S_inv;     % 反转
        k_S1;      % x驱动上升系数
        k_S2;      % x驱动失活系数
        k_S3;      % 恢复系数

        % ===== F的动力学参数 =====
        F_rest;    % 静息态
        F_min;     % 彻底失活
        F_max;     % 饱和上限
        F_thr;     % F上升阈值
        F_inv;     % 反转
        k_F1;      % x驱动上升系数
        k_F2;      % x驱动下降系数
        k_F3;      % 恢复系数
    end

    methods
        %% 初始化
        function obj = wendling(E0, S0, F0, dt, if_dy)
            % 构造函数
            obj = obj@base();
            obj.name = 'Wendling + Dynamic ESF';
            obj.params_size = 16; % 参数量
            obj.state_size = 10;  % 原文 y0..y9：四群体 + 慢抑制子群 y4，共五组 (y,x)
            obj.input_size = 1; % 输入维度
            obj.output_size = 1;  % 输出维度
            obj.dt = dt;

            % 初始化
            obj.current = 0.001*rand(obj.state_size,1);
            obj.history = zeros(obj.state_size,1);
            obj.input = zeros(obj.input_size,1);
            obj.output = zeros(obj.output_size,1);

            % 是否使用动态变量
            if if_dy == 1
                obj.if_dy = 1;
                obj.dy_size = 3; % 3个动态参数E、S、F

                obj.dy_ESF = [E0; S0; F0];


                % ===== E的动力学参数（按你的测试结果调）=====
                obj.E_rest = obj.dy_rest(1);
                obj.E_min  = obj.dy_min(1);
                obj.E_max  = obj.dy_max(1);
                obj.E_thr  = obj.dy_thr(1);
                obj.E_inv  = obj.dy_inv(1);
                obj.k_E1   = obj.x_in(1);
                obj.k_E2   = obj.x_de(1);
                obj.k_E3   = obj.dy_re(1);

                % ===== S的动力学参数 =====
                obj.S_rest = obj.dy_rest(2);
                obj.S_min  = obj.dy_min(2);
                obj.S_max  = obj.dy_max(2);
                obj.S_thr  = obj.dy_thr(2);
                obj.S_inv  = obj.dy_inv(2);
                obj.k_S1   = obj.x_in(2);
                obj.k_S2   = obj.x_de(2);
                obj.k_S3   = obj.dy_re(2);

                % ===== F的动力学参数 =====
                obj.F_rest = obj.dy_rest(3);
                obj.F_min  = obj.dy_min(3);
                obj.F_max  = obj.dy_max(3);
                obj.F_thr  = obj.dy_thr(3);
                obj.F_inv  = obj.dy_inv(3);
                obj.k_F1   = obj.x_in(3);
                obj.k_F2   = obj.x_de(3);
                obj.k_F3   = obj.dy_re(3);
            end

            obj.G_e = E0;  % E初始值
            obj.G_s = S0;  % S初始值
            obj.G_f = F0;  % F初始值

        end

        %% 关键微分过程
        function [d_ESF, dy] = derivatives(obj, y, input)
            % ---------- 1. 提取状态 ----------
            % 提取原有状态
            y0 = y(1); y1 = y(2);
            y2 = y(3); y3 = y(4);
            y4 = y(5); y5 = y(6);
            y6 = y(7); y7 = y(8);
            y8 = y(9); y9 = y(10);

            % S 型函数（本步要用多次，先建句柄）
            Sfun = @(v) 2*obj.e0./(1+exp(obj.r*(obj.v0-v)));

            % 模型输出：三个中间神经元群体 PSP 的净和，也是增益动力学的驱动来源
            y_out = y1 - y2 - y3;

            % ---------- 2. 增益的动力学 ----------
            % 注意：derivatives 只计算导数，绝不修改对象状态（修掉 E2：
            % 旧版在这里改 E_inv / k_E1，而 RK4 一个时间步会调用本函数 4 次）。
            if obj.if_dy == 1

                Gv = obj.dy_ESF(:);            % 当前 [G_e; G_s; G_f] (mV)

                if strcmpi(obj.gain_law, 'bistable')
                    % ===== 统一双稳形式（见 gain_bistable.m）=====
                    % 驱动 h = 归一化发放率，不整流、不减阈值（修掉 E1）
                    h     = Sfun(y_out) / (2*obj.e0);
                    d_ESF = gain_bistable(Gv, h, obj.bs);
                    Ge = Gv(1); Gs = Gv(2); Gf = Gv(3);
                else
                    % ===== legacy：分段线性 + 三个独立状态机（保留以复现旧图）=====
                    E = Gv(1);
                    S = Gv(2);
                    F = Gv(3);

                    % 计算有效刺激x, 输出定义：y_out = y0
                    % 背景放电基线（mV）：真实皮层不静默，静息 PSP 也会偶有阈上成分。
                    % 0 = 旧行为（静息 x 恒为 0，增益一动不动）；>0 时静息 x 恒正，
                    % G_e 会缓慢爬向平衡点 rest + k1·x0/k3 = 3.25 + 200·x0，
                    % 长时程里因此能出现"发作—恢复—再发作"。逐节点可取不同值。
                    x = max(0, y0 - obj.y_thr) + obj.gain_bias;

                    % E：正反馈，受x增加，本身有回落趋势
                    if E > obj.E_thr
                        obj.E_inv = 1;
                    end
                    if obj.E_inv == 0
                        dE =   obj.k_E1 * x - obj.k_E3 * (E - obj.E_rest);
                        dE = min(0.1,dE);
                    else
                        dE = - obj.k_E2 * x + obj.k_E3 * (obj.E_rest - E);
                    end

                    % S：先升后降，E超阈值后自动翻转
                    if S > obj.S_thr
                        obj.S_inv = 1;
                    end
                    if obj.S_inv == 0
                        dS =   obj.k_S1 * x - obj.k_S3 * (S - obj.S_rest);
                    else
                        dS = - obj.k_S2 * x + obj.k_S3 * (obj.S_rest - S);
                    end

                    % F：受S抑制，S降后延迟下跌
                    if F > obj.F_thr
                        obj.F_inv = 1;
                    end
                    if obj.F_inv == 0
                        dF =   obj.k_F1 * x - obj.k_F3 * (F - obj.F_rest);
                    else
                        dF = - obj.k_F2 * x + obj.k_F3 * (obj.F_rest - F);
                    end

                    % 失活后开始恢复，x驱动失效
                    if S == obj.S_min
                        obj.S_inv = 0;obj.k_S1 = 0;
                    end
                    if E == obj.E_min
                        obj.E_inv = 0;obj.k_E1 = 0;
                    end
                    if F == obj.F_min
                        obj.F_inv = 0;obj.k_F1 = 0;
                    end

                    % 达到静息态恢复x驱动
                    % gain_reset_tol = 0 时是旧行为（精确相等）；>0 时用容差，
                    % 使长时程仿真能出现下一次发作（见属性区注释）。
                    if obj.gain_reset_tol > 0
                        back_to_rest = abs(S - obj.S_rest) < obj.gain_reset_tol && ...
                                       abs(E - obj.E_rest) < obj.gain_reset_tol && ...
                                       abs(F - obj.F_rest) < obj.gain_reset_tol;
                    else
                        back_to_rest = (S == obj.S_rest && E == obj.E_rest && F == obj.F_rest);
                    end
                    if back_to_rest
                        obj.k_E1 = obj.x_in(1);obj.k_S1 = obj.x_in(2);obj.k_F1 = obj.x_in(3);
                    end

                    d_ESF = [dE; dS; dF];
                    Ge = E; Gs = S; Gf = F;
                end

            % 不使用动态参数
            else
                d_ESF = zeros(3,1);
                Ge = obj.G_e; Gs = obj.G_s; Gf = obj.G_f;
            end

            % 计算输入到各群体（每个 v 就是原文 S 括号里的那一串）
            v1 = y1 - y2 - y3;                  % PC  : y1 - y2 - y3       （无 C，原文如此）
            v2 = obj.C_ep*y0;                   % EIN : C1 y0
            v3 = obj.C_sp*y0;                   % SIN : C3 y0
            v4 = obj.C_fp*y0 - obj.C_pf*y4;     % FIN : C5 y0 - C6 y4     （C6 用的是 y4 不是 y3）
            v5 = obj.C_sp*y0;                   % SIN': C3 y0              （与 v3 完全相同）

            % 计算脉冲密度
            z1 = Sfun(v1);
            z2 = Sfun(v2);
            z3 = Sfun(v3);
            z4 = Sfun(v4);
            z5 = Sfun(v5);

            % 微分方程（逐行对照 Wendling 2002 式(1)，A→G_e，B→G_s，G→G_f，a→w_e，b→w_s，g→w_f）
            %   dy(1)  = y5                          y0' = y5
            %   dy(6)  = G_e w_e  z1  - 2 w_e y5 - w_e^2 y0     ← A a S(y1-y2-y3)
            %   dy(2)  = y6                          y1' = y6
            %   dy(7)  = G_e w_e (input + C_pe z2) ...          ← A a {p(t) + C2 S(C1 y0)}
            %   dy(3)  = y7                          y2' = y7
            %   dy(8)  = G_s w_s C_ps z3 ...                    ← B b C4 S(C3 y0)
            %   dy(4)  = y8                          y3' = y8
            %   dy(9)  = G_f w_f C_fs z4 ...                    ← G g C7 S(C5 y0 - C6 y4)
            %   dy(5)  = y9                          y4' = y9
            %   dy(10) = G_s w_s z5 ...                         ← B b S(C3 y0)

            % 初始化微分方程
            dy = zeros(obj.state_size,1);

            dy(1) = y5;
            dy(6) = Ge*obj.w_e*z1 - 2*obj.w_e*y5 - obj.w_e^2*y0;

            dy(2) = y6;
            dy(7) = Ge*obj.w_e*(input+obj.C_pe*z2) - 2*obj.w_e*y6 - obj.w_e^2*y1;

            dy(3) = y7;
            dy(8) = Gs*obj.w_s*obj.C_ps*z3 - 2*obj.w_s*y7 - obj.w_s^2*y2;

            dy(4) = y8;
            dy(9) = Gf*obj.w_f*obj.C_fs*z4 - 2*obj.w_f*y8 - obj.w_f^2*y3;

            dy(5) = y9;
            dy(10) = Gs*obj.w_s*z5 - 2*obj.w_s*y9 - obj.w_s^2*y4;

        end

        %% 积分
        function [dy_ESF, y, output] = dynamics(obj, y, input)
            % 四阶Runge-Kutta算法求解一步
            % input: 外部输入（高斯白噪声）

            [l1,k1] = obj.derivatives(y, input);
            [l2,k2] = obj.derivatives(y + 0.5 * obj.dt * k1, input);
            [l3,k3] = obj.derivatives(y + 0.5 * obj.dt * k2, input);
            [l4,k4] = obj.derivatives(y + obj.dt * k3, input);

            % 更新状态
            if obj.if_dy==1
                obj.dy_ESF = obj.dy_ESF + (obj.dt / 6) * (l1 + 2*l2 + 2*l3 + l4);
                % 范围限制：dy_min/dy_max 就是三个增益的物理区间 [3,9]/[3,31]/[5,18]，
                % 双稳形式的 φ∈[0,1] 与之逐位对应，所以这里只是数值兜底，
                % 不是"钳位式"的建模假设（旧版的钳位是在 derivatives 里按 if 做的）。
                obj.dy_ESF = max(obj.dy_min, min(obj.dy_ESF, obj.dy_max));
                % 同步三个增益属性，方便外部查看/画图（derivatives 不再改它们）
                obj.G_e = obj.dy_ESF(1);
                obj.G_s = obj.dy_ESF(2);
                obj.G_f = obj.dy_ESF(3);
            end
            y = y + (obj.dt / 6) * (k1 + 2*k2 + 2*k3 + k4);

            % 范围限制
            % obj.dy_ESF = max(obj.dy_min, min(obj.dy_ESF, obj.dy_max));
            % obj.dy_ESF(1) = max(obj.E_min, min(obj.dy_ESF(1), obj.E_max));
            % obj.dy_ESF(2) = max(obj.S_min, min(obj.dy_ESF(2), obj.S_max));
            % obj.dy_ESF(3) = max(obj.F_min, min(obj.dy_ESF(3), obj.F_max));

            % 计算输出（锥体细胞的突触后电位）
            %output = y(1);
            output = y(2) - y(3) - y(4);
            dy_ESF = obj.dy_ESF;
        end

        %% 更改参数

        function obj = ESF(obj, E_t, S_t, F_t)
            % ESF  改写三个突触增益。
            %
            %   静态模式（if_dy = 0）下增益是常数，只能写 G_e/G_s/G_f；
            %   旧版无条件写 dy_ESF，而静态模型的 dy_size = 0，于是下一个整合步
            %   System 里 `dy_param(idx) = dy_ESF` 会因左右元素数不同而报错
            %   （只能在 t > esf_time 之后才炸，所以 5 s 能跑通、20 s 会崩）。
            %   这是做 (B,G) 静态参数扫描时必踩的坑，这里按模式分流。
            if obj.if_dy == 1
                obj.dy_ESF = [E_t; S_t; F_t];
            else
                obj.G_e = E_t; obj.G_s = S_t; obj.G_f = F_t;
            end
        end

        function obj = set_gain_law(obj, law, bs)
            % SET_GAIN_LAW 切换增益动力学形式（'legacy' / 'bistable'）
            %
            %   obj = obj.set_gain_law('bistable')
            %   obj = obj.set_gain_law('bistable', bs_struct)
            %
            % 由 simulate_network 从 cfg.gain_law / cfg.bs_params 调用；
            % bs 省略时沿用属性里的默认参数（见 properties 段的 bs）。
            if nargin >= 2 && ~isempty(law)
                obj.gain_law = lower(law);
            end
            if nargin >= 3 && ~isempty(bs)
                f = fieldnames(bs);
                for k = 1:numel(f)
                    if isfield(obj.bs, f{k})
                        obj.bs.(f{k}) = bs.(f{k})(:);
                    else
                        error('未知的双稳参数名：%s', f{k});
                    end
                end
            end
        end

        function obj = put(obj, params)
            % 提取参数
            if isfield(params, 'G_e'), obj.G_e = params.G_e; end
            if isfield(params, 'G_s'), obj.G_s = params.G_s; end
            if isfield(params, 'G_f'), obj.G_f = params.G_f; end
            if isfield(params, 'w_e'), obj.w_e = params.w_e; end
            if isfield(params, 'w_s'), obj.w_s = params.w_s; end
            if isfield(params, 'w_f'), obj.w_f = params.w_f; end
            if isfield(params, 'C_pe'), obj.C_pe = params.C_pe; end
            if isfield(params, 'C_ps'), obj.C_ps = params.C_ps; end
            if isfield(params, 'C_pf'), obj.C_pf = params.C_pf; end
            if isfield(params, 'C_ep'), obj.C_ep = params.C_ep; end
            if isfield(params, 'C_sp'), obj.C_sp = params.C_sp; end
            if isfield(params, 'C_fp'), obj.C_fp = params.C_fp; end
            if isfield(params, 'C_fs'), obj.C_fs = params.C_fs; end
            if isfield(params, 'e0'), obj.e0 = params.e0; end
            if isfield(params, 'r'), obj.r = params.r; end
            if isfield(params, 'v0'), obj.v0 = params.v0; end
        end

    end
end