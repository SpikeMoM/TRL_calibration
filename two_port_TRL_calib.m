clear;
clc;
close all;
%% 输入端口
Prot_in = input('\nNumber of Input Port: Prot_in = ');
Prot_out = input('\nNumber of Output Port: Prot_out = ');
%% 计算误差系数
point = 201;
%%导入TRL校准数据
S_T = dlmread(['S_Thru_reim.txt'],'',[7 0 point+6 4]);
S_R = dlmread(['S_Refl_reim.txt'],'',[7 0 point+6 2]);
S_L = dlmread(['S_Line_reim.txt'],'',[7 0 point+6 4]);
%%数据处理
freq = S_T(:, 1);
R11 = complex( S_R(:, 2), S_R(:, 3));
L11 = complex( S_L(:, 2), S_L(:, 3));
L12 = complex( S_L(:, 4), S_L(:, 5));
T11 = complex( S_T(:, 2), S_T(:, 3));
T12 = complex( S_T(:, 4), S_T(:, 5));
%%计算
line_delay = (L12.^2 + T12.^2 - (T11 - L11).^2 + ...
    power(((L12.^2 + T12.^2 - (T11 - L11).^2).^2) - 4*L12.^2.*T12.^2, 0.5))./(2*L12.*T12);
line_delay_real = real(line_delay);
line_delay_imag = imag(line_delay);
line_delay = complex(line_delay_real, -abs(line_delay_imag));
line_angle = angle(line_delay)*180/pi;
S22 = (T11 - L11)./(T12 - L12.*line_delay);
S11 = T11 - S22.*T12;
S12 = power((T12.*(1-S22.^2)), 0.5);
RL =(R11 - S11)./(S12.^2 + S22.*(R11 - S11)); 
[eA, eB, eC, eD] = S2ABCD(S11, S12, S12, S22);
eABCD = zeros(2,2,point);
inv_eABCD = zeros(2,2,point);
eABCD(1, 1, :) = eA;
eABCD(1, 2, :) = eB;
eABCD(2, 1, :) = eC;
eABCD(2, 2, :) = eD;
for i = 1:point
    inv_eABCD(:, :, i) = inv(eABCD(:, :, i));
end
%% 校准
%%导入DUT测试数据
S_eDUT = dlmread(['DUT_meas_reim.txt'],'',[7 0 point+6 4]);
%%数据处理
S_eDUT_11 = complex( S_eDUT(:, 2), S_eDUT(:, 3));
S_eDUT_21 = complex( S_eDUT(:, 4), S_eDUT(:, 5));
[eDUT_A, eDUT_B, eDUT_C, eDUT_D] = ...
    S2ABCD(S_eDUT_11, S_eDUT_21, S_eDUT_21, S_eDUT_11);
eDUT_ABCD = zeros(2,2,point);
DUT_ABCD_calib = zeros(2,2,point);
eDUT_ABCD(1, 1, :) = eDUT_A; eDUT_ABCD(1, 2, :) = eDUT_B;
eDUT_ABCD(2, 1, :) = eDUT_C; eDUT_ABCD(2, 2, :) = eDUT_D;
%%级联误差和测试ABCD矩阵得到器件的ABCD矩阵
for i = 1:point
    DUT_ABCD_calib(:, :, i) = inv_eABCD(:, :, i)*eDUT_ABCD(:, :, i)*eABCD(:, :, i);
end
DUT_S_calib = zeros(2,2,point);
%%将器件的ABCD矩阵转化为S参数矩阵
for i = 1:point
    [DUT_S_calib(1, 1, i), DUT_S_calib(1, 2, i), DUT_S_calib(2, 1, i), DUT_S_calib(2, 2, i)] = ...
        ABCD2S(DUT_ABCD_calib(1, 1, i),DUT_ABCD_calib(1, 2, i),DUT_ABCD_calib(2, 1, i),DUT_ABCD_calib(2, 2, i));
end
DUT_S11_calib = reshape(DUT_S_calib(1,1,:),point,1);
DUT_S12_calib = reshape(DUT_S_calib(1,2,:),point,1);
DUT_S21_calib = reshape(DUT_S_calib(2,1,:),point,1);
DUT_S22_calib = reshape(DUT_S_calib(2,2,:),point,1);%得到器件平面的计算结果
%% 比较
f_low = 0;
f_up = 24;
f_step = 2;
y_low = -60;
y_up = 0;
y_step = 10;
%%导入仿真数据
S_sim = dlmread(['DUT_sim_dB.txt'],'',[7 0 point+6 2]);
DUT_S11_sim = S_sim(:, 2);
DUT_S21_sim = S_sim(:, 3);
export_S11 = 20*log10( abs(DUT_S11_calib));
export_S21 = 20*log10( abs(DUT_S21_calib));
%% 绘图
figure('Position', [100, 100, 600, 400]);
b1 = plot(freq, DUT_S11_sim, 'r-', freq, 20*log10( abs(DUT_S11_calib)), 'm--',...
    freq, DUT_S21_sim, 'b-', freq, 20*log10( abs(DUT_S21_calib)), 'c--');
set( gca, 'linewidth', 1);
set( b1, 'linewidth', 2);
axis([f_low f_up y_low y_up]);
xlabel('Frequency(GHz)');
ylabel('S (dB)');
legend(['S_{',num2str(Prot_in),num2str(Prot_in),'}(sim)'], ...
    ['S_{',num2str(Prot_in),num2str(Prot_in),'}(cal)'], ...
    ['S_{',num2str(Prot_out),num2str(Prot_in),'}(sim)'], ...
    ['S_{',num2str(Prot_out),num2str(Prot_in),'}(cal)'], 'location', 'southeast');
grid off;
ax = gca;
ax.FontSize = 16;
ax.FontName = 'times new roman';
ax.XTick = f_low:f_step:f_up;
ax.YTick = y_low:y_step:y_up;






