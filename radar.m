%===========================================================================================
%????????????该程序完成16个脉冲信号的【脉压、动目标显示/动目标检测（MTI/MTD）】?????
%===========================================================================================
%??程序中根据每个学生学号的末尾三位（依次为XYZ）来决定仿真参数，034?
%??目标距离为[3000?8025?9000+(Y*10 + Z)*200?8025]，4个目标
%??目标速度为[50?0?(Y*10 + X + Z)*6?100]
%?===================================================================================%?
%????????????????????????????????????雷达参数???????????????????????????????????????%?
%?===================================================================================%?
clear
close all
C = 3.0e8;%光速(m/s)?
RF = 3.140e9/2;%雷达射频?1.57GHz?
Lambda = C/RF;%雷达工作波长?
PulseNumber = 16;%回波脉冲数?
BandWidth = 2.0e6;%发射信号带宽?带宽B = 1/τ，τ是脉冲宽度
PulseWidth = 42.0e-6;%发射信号时宽?
PRT = 240e-6;%?雷达发射脉冲重复周期(s),240us对应1/2*240*300 = 36000米最大无模糊距离
PRF = 1/PRT;%雷达发射脉冲频率
Fs = 2.0e6;%采样频率?
NoisePower = -12;%(dB);%噪声功率（目标为0dB）?
%?---------------------------------------------------------------
SampleNumber = fix(Fs*PRT);%计算一个脉冲周期的采样点数480 = 240 * 2 
TotalNumber = SampleNumber*PulseNumber;%总的采样点数480*16 =
BlindNumber = fix(Fs*PulseWidth);%计算一个脉冲周期的盲区-遮挡样点数
%===================================================================================%
%????????????????????????????????????目标参数????????????????????????????????????????
%===================================================================================%

TargetNumber = 1;%目标个数
SigPower(1:TargetNumber) = [1 1 1 0.25];%目标功率,无量纲
TargetDistance(1:TargetNumber) = [3000 8025 15800 8025];%目标距离,单位m???距离参数为[3000?8025?9000+(Y*10 + Z)*200?8025]
DelayNumber(1:TargetNumber) = fix(Fs*2*TargetDistance(1:TargetNumber)/C);%?把目标距离换算成采样点（距离门）?fix函数向0靠拢取整
TargetVelocity(1:TargetNumber) = [50 0 204 100];%目标径向速度?单位m/s???速度参数为[50?0?(Y*10 + X + Z)*6?100]
TargetFd(1:TargetNumber) = 2*TargetVelocity(1:TargetNumber)/Lambda;%计算目标多卜勒频移2v/λ


%====================================================================================%
%???????????????????????????????????产生线性调频信号??????????????????????????????????%
%====================================================================================%
number = fix(Fs*PulseWidth);%回波的采样点数 = 脉压系数长度 = 暂态点数目 + 1

if rem(number,2)~=0 %rem求余
    number = number + 1;
end %把number变为偶数

for i = - fix(number/2):fix(number/2) - 1
    Chirp(i + fix(number/2)+1) = exp(j*(pi*(BandWidth/PulseWidth)*(i/Fs)^2));产生复数矩阵Chirp
end

coeff = conj(fliplr(Chirp));%把Chirp矩阵翻转并把复数共轭，产生脉压系数
figure(1);%脉压系数的实部
plot(real(Chirp));
axis([0 90 -1.5 1.5]);
title('脉压系数实部');

%-----------------------------产生目标回波串-----------------%
SignalAll = zeros(1,TotalNumber);%所有脉冲的信号,先填0

for k = 1:TargetNumber %依次产生各个目标
    SignalTemp = zeros(1,SampleNumber);%一个PRT
    SignalTemp(DelayNumber(k) + 1:DelayNumber(k) + number) = sqrt(SigPower(k))*Chirp;

    %一个脉冲的1个目标（未加多普勒速度）(DelayNumber(k)+1):(DelayNumber(k)+number)
    Signal = zeros(1,TotalNumber);

    for i = 1:PulseNumber%?16个回波脉冲
       Signal((i-1)*SampleNumber + 1:i*SampleNumber) = SignalTemp;%每个目标把16个SignalTemp排在一起
    end

    FreqMove = exp(j*2*pi*TargetFd(k)*(0:TotalNumber - 1)/Fs);%目标的多普勒速度*时间 = 目标的多普勒相移
    Signal = Signal.*FreqMove;%加上多普勒速度后的16个脉冲1个目标
    SignalAll = SignalAll + Signal;%加上多普勒速度后的16个脉冲4个目标
end

%-------------------------产生第4个目标的回波串-------%
%fi = pi/3;
%SignalTemp = zeros(1,SampleNumber);%一个脉冲
%SignalTemp(DelayNumber(4) + 1:DelayNumber(4) + number) = sqrt(SigPower(4))*exp(j*fi)*Chirp;%一个脉冲的1个目标(未加多普勒速度)
%Signal = zeros(1,TotalNumber);
%
%for i = 1:PulseNumber
%    Signal((i-1)*SampleNumber + 1:i*SampleNumber) = SignalTemp;
%end
%
%FreqMove = exp(j*2*pi*TargetFd(4)*(0:TotalNumber-1)/Fs);%目标的多普勒速度*时间 = 目标的多普勒相移
%Signal = Signal.*FreqMove;
%SignalAll = SignalAll + Signal;



figure(2);
subplot(2,1,1);plot(real(SignalAll),'r-');
title('目标信号的实部');
grid on;
zoom on;
subplot(2,1,2);plot(imag(SignalAll));
title('目标信号的虚部');
grid on;
%zoom on;
%====================================================================================%?
%???????????????????????????????????产生系统噪声信号??????????????????????????????????%?
%====================================================================================%?
SystemNoise = normrnd(0,10^(NoisePower/10),1,TotalNumber)+j*normrnd(0,10^(NoisePower/10),1,TotalNumber);
%均值为0，标准差为10^(NoisePower/10)的噪声?
%====================================================================================%?
%???????????????????????????????????总的回波信号?????????????????????????????????????%?
%====================================================================================%?
Echo = SignalAll + SystemNoise;%?+SeaClutter + TerraClutter，加噪声之后的回波?
for i = 1:PulseNumber%在接收机闭锁期,接收的回波为0?
    Echo((i-1)*SampleNumber + 1:(i-1)*SampleNumber + number)=0;%发射时接收为0?
end
figure(3);%加噪声之后的总回波信号?
subplot(2,1,1);plot(real(Echo),'r-');title('总回波信号的实部,闭锁期为0');
subplot(2,1,2);plot(imag(Echo));title('总回波信号的虚部,闭锁期为0');



%================================时域脉压=================================%
pc_time0 = conv(Echo,coeff);%pc_time0为Echo和coeff的卷积?
pc_time1 = pc_time0(number:TotalNumber + number-1);%去掉暂态点?number-1个?
figure(4);%时域脉压结果的幅度?
subplot(2,1,1);plot(abs(pc_time0),'r-');title('时域脉压结果的幅度,有暂态点');%pc_time0的模的曲线?
subplot(2,1,2);plot(abs(pc_time1));title('时域脉压结果的幅度,无暂态点');%pc_time1的模的曲线?




%?================================频域脉压=================================%?
Echo_fft = fft(Echo,8192);%理应进行TotalNumber + number-1点FFT,但为了提高运算速度,进行了8192点的FFT?
coeff_fft = fft(coeff,8192);
pc_fft = Echo_fft.*coeff_fft;
pc_freq0 = ifft(pc_fft);
figure(5);
subplot(2,1,1);plot(abs(pc_freq0(1:TotalNumber + number-1)));
title('频域脉压结果的幅度,有前暂态点');
subplot(2,1,2);plot(abs(pc_time0(1:TotalNumber + number-1)-pc_freq0(1:TotalNumber + number-1)),'r');
title('时域和频域脉压的差别');
pc_freq1 = pc_freq0(number:TotalNumber + number-1);%去掉暂态点?number-1个,后填充点若干(8192-number + 1-TotalNumber)?

%?================按照脉冲号、距离门号重排数据=================================%?
for i = 1:PulseNumber
    pc(i,1:SampleNumber)=pc_freq1((i-1)*SampleNumber + 1:i*SampleNumber);%每个PRT为一行，每行480个采样点的数据?
end
figure(6);
plot(abs(pc(1,:)));
title('频域脉压结果的幅度,没有暂态点');
%?================MTI（动目标显示）,对消静止目标和低速目标---可抑制杂波=================================%?
for i = 1:PulseNumber-1 %滑动对消，少了一个脉冲????
    mti(i,:)=pc(i + 1,:)-pc(i,:);
end
figure(7);
mesh(abs(mti));
title('MTI??result');
%?================MTD（动目标检测）,区分不同速度的目标，有测速作用=================================%?
mtd = zeros(PulseNumber,SampleNumber);
for i = 1:SampleNumber
    buff(1:PulseNumber)=pc(1:PulseNumber,i);
    buff_fft = fft(buff);
    mtd(1:PulseNumber,i)=buff_fft(1:PulseNumber);
end
figure(8);
mesh(abs(mtd));
title('MTD??result');
%=======================================虚实矩阵转换========================================%
coeff_fft_c = zeros(1,2*8192);
for i = 1:8192
    coeff_fft_c(2*i-1)=real(coeff_fft(i));
    coeff_fft_c(2*i)=imag(coeff_fft(i));
end
echo_c = zeros(1,2*TotalNumber);
for i = 1:TotalNumber
    echo_c(2*i-1)=real(Echo(i));
    echo_c(2*i)=imag(Echo(i));
end

