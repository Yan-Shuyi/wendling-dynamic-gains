classdef base < handle
% BASE_MODEL  Abstract base class for neural population models.
%
%   classdef base < handle
%
% 作用
%   给所有神经群体模型定一个统一接口，System.m 只认这个接口，不关心具体是
%   哪个模型。想再加一种群体模型（比如 Jansen-Rit），只要继承 base 并实现
%   simulate / derivatives 即可，网络层一行都不用改。
%
% 子类必须提供的属性
%   name         模型名（字符串）
%   params_size  参数个数
%   state_size   状态变量个数（wendling 为 10：5 个群体 × (y, x)）
%   input_size   输入维度
%   output_size  输出维度（wendling 为 1，即 y_out）
%
% 子类必须实现的方法
%   simulate(input)  推进一个时间步，返回 [dy_ESF, x, y]
%
% 调用关系
%   下游：wendling（唯一的实现）
%   上游：System.network_dynamics
    properties
        name            % 模型名称
        params_size     % 模型参数量
        state_size      % 状态变量维度
        current         % 当前状态变量
        history         % 历史状态记录
        input_size      % 输入维度
        input           % 当前输入
        output_size     % 输出维度
        output          % 当前输出
    end
    
    methods (Abstract)
        [dy_ESF, x, output] = dynamics(obj, x, input)
        % DYNAMICS 定义模型微分方程 (抽象方法)
        % 输入:
        %   x - 当前状态向量
        %   input - 外部输入向量
        % 输出:
        %   x - 状态
    end
    
    methods
        function obj = base()
        end
        
        function [dy_ESF, x, output] = simulate(obj, input)
            % SIMULATE 执行模型仿真
            % 输入:
            %   x - 当前状态
            %   inputs - 外部输入函数句柄 inputs
            
            obj.input = input;

            % 封装动力学方程
            [dy_ESF, x, output] = dynamics(obj, obj.current, input);

            % 记录历史状态
            obj.history = obj.current;

            % 记录当前状态
            obj.current = x;
            obj.output = output;

        end
    end
end