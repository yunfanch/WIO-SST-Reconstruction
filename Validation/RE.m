function RE_index = RE( x1, x2, cal, ver)
% The RE statistic tests whether a reconstruction provides a better estimate 
% of climatic variability than climatology (i.e., the mean of the meteorological 
% data in the calibration period [Cook and Kairiukstis, 1990; Cook et al., 1994]).
% x1 is the actual data
% x2 is the estimated data
% cal is the calibration period
% ver is the verification period
    ver_sum = sum((x1(ver)-x2(ver)).^2,'omitmissing');
    ver_sum2 = sum((x1(ver)-mean(x1(cal),'omitmissing')).^2,'omitmissing');
    RE_index = 1 - ver_sum./ver_sum2;
end

