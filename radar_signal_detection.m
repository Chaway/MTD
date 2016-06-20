% This M-file shows the radar signal detection including MIT, MTD and CFAR.

clear all;
close all;
clc;

%------------------------------��������------------------------------------------------------------------------------------------------------
fr = 200;           % PRI��Hz��
tr = 1 / fr;        % ��ǰ���ڣ��������ظ����ڣ�s��
fd1 = 90;           % ��Ŀ��1�Ķ�����Ƶ�ƣ�Hz��
fd2 = 140;          % ��Ŀ��2�Ķ�����Ƶ�ƣ�Hz��
M = 6700;           % �ز��źŵľ������
N = 10;             % ��ǰ������
n_reference = 11;   % CFAR���쵥Ԫÿ��Ĳο���Ԫ��
n_protect = 2;      % CFAR���쵥Ԫÿ��ı�����Ԫ��
alpha = 2.9;        % CFAR��������

%-----------------------------����N����ǰ���ڣ�ÿ����ǰ���ں�M�������ĸ����ݣ�����ֻ������Ŀ�꣬���Ӳ�------------------------------------------
for m = 1 : N
    moving_target1(m) = 10 * exp(j * 2 * m * pi * fd1 * tr);    % ��m����ǰ���ڣ���Ŀ��1��ĳ������ڵĸ�����
    moving_target2(m) = 10 * exp(j * 2 * m * pi * fd2 * tr);    % ��m����ǰ���ڣ���Ŀ��2��ĳ������ڵĸ�����
end

x = ones(1, 3);     % ÿ����Ŀ��ռ�����������
n = length(x);

for m = 1 : N
    sequence(m, :) = [[zeros(1, 1700), x, zeros(1, 3400 - 1700 - n)] * moving_target1(m), [zeros(1, 1600), x, zeros(1, M - 3400 - 1600 - n)] * moving_target2(m)];
                                                                   % ��m����ǰ���ڣ��Ŀ��1ռ��1701��1702��1703��������㣬�Ŀ��2ռ�ݵ�5001��5002��5003���������
end

%------------------------------�����Ӳ�����������뵽�������ڵ��������У��γɻز�����ź�--------------------------------------------------------
for m = 1 : N
    clutter_i(m, :) = randn(1, M);
    clutter_q(m, :) = randn(1, M);
    clutter(m, :) = 0.5 * (clutter_i(m, :) + j * clutter_q(m, :));   % �Ӳ�
end

for m = 1 : N
    mix_signal(m, :) = sequence(m, :) + clutter(m, :);        % ���ɻ���ź�
end

%------------------------------�����������MTI��----------------------------------------------------------------------------------------------
for m = 1 : N - 2
    mti(m, :) = mix_signal(m, :) - 2 * mix_signal(m + 1, :) + mix_signal(m + 2, :);     % �������������
end

%------------------------------FFT���㣨MTD��------------------------------------------------------------------------------------------------
for n = 1 : M
    mtd(:, n) = fft(mti(:, n) .* hamming(N - 2, 'periodic'), N - 2);       % N-2��FFT����
end

mtd_amplitude = abs(mtd);       % MTD��ͨ��ֱ�����ģֵ

set(0, 'defaultfigurecolor', 'w');
for m = 1 : N - 2
    figure(m + 5);
    plot(mtd_amplitude(m, :), 'r');
    grid on;
    xlim([0, M]);
    xlabel('�����');
    title(['MTD�� ' num2str(m - 1) ' ͨ��������']);
end

%-----------GOCA CFAR---------------------------------------------------------------------------------------------------
for m = 1 : N - 2                     % �ֱ��MTD N-2��ͨ���ڵľ�����������޼��
    
    z = mtd_amplitude(m, :);
    
    for n = 1 : 1 + n_protect         % ���쵥Ԫ
        b2(n) = mean(z(n + 1 + n_protect : n + n_protect + n_reference));  % ���쵥Ԫ�Ҳ�Ĳο������ֵ
        d(n) = b2(n);           
    end
    
    for n = 1 + n_protect + 1 : n_reference + n_protect + 1
        b1(n) = mean(z(1 : n - 1 - n_protect));                            % ���쵥Ԫ���Ĳο������ֵ
        b2(n) = mean(z(n + 1 + n_protect : n + n_protect + n_reference));
        d(n) = max(b1(n), b2(n));                                          % ����ο�����ֵѡ��
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
        
        T(n) = alpha * d(n);                % ��n�������ļ������
        
        if z(n) > T(n)                      % ������
            s_cfar(n) = z(n) - T(n);        % ���þ�����������Ϊ��ֵ
        end
    end
    
    temp(m, :) = s_cfar;                    % �洢����ͨ���ļ����
    
    figure(m + 40);
    plot(z, 'r');
    hold on;
    plot(T, 'g');
    xlabel('�����');
    grid on;
    xlim([0, M]);
    title(['MTD�� ' num2str(m - 1) ' ͨ���������趨']);
    legend('����ź�', '����');
    
    figure(m + 80);
    hold off;
    plot(s_cfar, 'r');
    xlabel('�����');
    grid on ;
    xlim([0, M]);
    title(['MTD�� ' num2str(m - 1) ' ͨ�������']);
end

for  n = 1 : M
    cfar(:, n) = max(temp(:, n));           % ȡN-2��ͨ���е����ֵ��Ϊ�þ����ļ����
end

figure(115);
plot(cfar, 'r');
title('GOCA CFAR�����');
xlabel('�����'); 
xlim([0, M]);
grid on ;

        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    