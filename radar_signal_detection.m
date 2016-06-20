% This M-file shows the radar signal detection including MIT, MTD and CFAR.

clear all;
close all;
clc;

%------------------------------参数设置------------------------------------------------------------------------------------------------------
fr = 200;           % PRI（Hz）
tr = 1 / fr;        % 导前周期，即脉冲重复周期（s）
fd1 = 90;           % 动目标1的多普勒频移（Hz）
fd2 = 140;          % 动目标2的多普勒频移（Hz）
M = 6700;           % 回波信号的距离点数
N = 10;             % 导前周期数
n_reference = 11;   % CFAR待检单元每侧的参考单元数
n_protect = 2;      % CFAR待检单元每侧的保护单元数
alpha = 2.9;        % CFAR门限因子

%-----------------------------生成N个导前周期，每个导前周期含M个距离点的复数据，其中只包含动目标，无杂波------------------------------------------
for m = 1 : N
    moving_target1(m) = 10 * exp(j * 2 * m * pi * fd1 * tr);    % 第m个导前周期，动目标1在某距离点内的复数据
    moving_target2(m) = 10 * exp(j * 2 * m * pi * fd2 * tr);    % 第m个导前周期，动目标2在某距离点内的复数据
end

x = ones(1, 3);     % 每个动目标占据三个距离点
n = length(x);

for m = 1 : N
    sequence(m, :) = [[zeros(1, 1700), x, zeros(1, 3400 - 1700 - n)] * moving_target1(m), [zeros(1, 1600), x, zeros(1, M - 3400 - 1600 - n)] * moving_target2(m)];
                                                                   % 第m个导前周期，令动目标1占据1701、1702和1703三个距离点，令动目标2占据第5001、5002和5003三个距离点
end

%------------------------------生成杂波，并将其加入到各个周期的数据序列，形成回波混合信号--------------------------------------------------------
for m = 1 : N
    clutter_i(m, :) = randn(1, M);
    clutter_q(m, :) = randn(1, M);
    clutter(m, :) = 0.5 * (clutter_i(m, :) + j * clutter_q(m, :));   % 杂波
end

for m = 1 : N
    mix_signal(m, :) = sequence(m, :) + clutter(m, :);        % 生成混合信号
end

%------------------------------三脉冲对消（MTI）----------------------------------------------------------------------------------------------
for m = 1 : N - 2
    mti(m, :) = mix_signal(m, :) - 2 * mix_signal(m + 1, :) + mix_signal(m + 2, :);     % 三脉冲对消运算
end

%------------------------------FFT运算（MTD）------------------------------------------------------------------------------------------------
for n = 1 : M
    mtd(:, n) = fft(mti(:, n) .* hamming(N - 2, 'periodic'), N - 2);       % N-2点FFT运算
end

mtd_amplitude = abs(mtd);       % MTD各通道直接输出模值

set(0, 'defaultfigurecolor', 'w');
for m = 1 : N - 2
    figure(m + 5);
    plot(mtd_amplitude(m, :), 'r');
    grid on;
    xlim([0, M]);
    xlabel('距离点');
    title(['MTD第 ' num2str(m - 1) ' 通道输出结果']);
end

%-----------GOCA CFAR---------------------------------------------------------------------------------------------------
for m = 1 : N - 2                     % 分别对MTD N-2个通道内的距离点做过门限检测
    
    z = mtd_amplitude(m, :);
    
    for n = 1 : 1 + n_protect         % 待检单元
        b2(n) = mean(z(n + 1 + n_protect : n + n_protect + n_reference));  % 待检单元右侧的参考窗求均值
        d(n) = b2(n);           
    end
    
    for n = 1 + n_protect + 1 : n_reference + n_protect + 1
        b1(n) = mean(z(1 : n - 1 - n_protect));                            % 待检单元左侧的参考窗求均值
        b2(n) = mean(z(n + 1 + n_protect : n + n_protect + n_reference));
        d(n) = max(b1(n), b2(n));                                          % 两侧参考窗均值选大
    end
    
    for n = 1 + n_reference + n_protect + 1 : M - n_reference - n_protect
        b1(n) = mean(z(n - n_protect - n_reference : n - 1 - n_protect));
        b2(n) = mean(z(n + 1 + n_protect : M));
        d(n) = max(b1(n), b2(n));
    end
    
    for n = M - n_reference - n_protect + 1 : M
        b1(n) = mean(z(n - n_protect - n_reference : n - 1 - n_protect));
        b2(n) = mean(z(n + 1 + n_protect : M));
        d(n) = max(b1(n), b2(n));
    end
    
    s_cfar = zeros(1, M);
    
    for n = 1 : M
        
        T(n) = alpha * d(n);                % 第n个距离点的检测门限
        
        if z(n) > T(n)                      % 过门限
            s_cfar(n) = z(n) - T(n);        % 将该距离点的数据置为差值
        end
    end
    
    temp(m, :) = s_cfar;                    % 存储各个通道的检测结果
    
    figure(m + 40);
    plot(z, 'r');
    hold on;
    plot(T, 'g');
    xlabel('距离点');
    grid on;
    xlim([0, M]);
    title(['MTD第 ' num2str(m - 1) ' 通道的门限设定']);
    legend('混合信号', '门限');
    
    figure(m + 80);
    hold off;
    plot(s_cfar, 'r');
    xlabel('距离点');
    grid on ;
    xlim([0, M]);
    title(['MTD第 ' num2str(m - 1) ' 通道检测结果']);
end

for  n = 1 : M
    cfar(:, n) = max(temp(:, n));           % 取N-2个通道中的最大值作为该距离点的检测结果
end

figure(115);
plot(cfar, 'r');
title('GOCA CFAR检测结果');
xlabel('距离点'); 
xlim([0, M]);
grid on ;

        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    