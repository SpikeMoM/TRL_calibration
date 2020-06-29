function [ S11, S12, S21, S22 ] = ABCD2S( A, B, C, D )
%ABCD矩阵转化为S参数
S11 = (A+B-C-D)./(A+B+C+D);
S12 = (2*A.*D-2*B.*C)./(A+B+C+D);
S21 = 2./(A+B+C+D);
S22 = (-A+B-C+D)./(A+B+C+D);
end

