function CE_index = CE( x1, x2, ver)
% The CE statistic is similar to RE except that its benchmark for determining 
% model skill is the verification period and not the calibration period (i.e., 
% the difference between RE and CE is in the denominator term).
% x1 is the actual data
% x2 is the estimated data
% cal is the calibration period
% ver is the verification period
    ver_sum = sum((x1(ver)-x2(ver)).^2,'omitmissing');
    cal_sum = sum((x1(ver)-mean(x1(ver),'omitmissing')).^2,'omitmissing');
    CE_index = 1 - ver_sum./cal_sum;
end

