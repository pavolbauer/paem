function [] = changetspan( path, tend )
%change_tspan Summary of this function goes here
%   Detailed explanation goes here

for i = [2 8:8:64]
    m = matfile (strcat(path, num2str(i), 'part.mat'), 'Writable', true);
    m.tspan = LINSPACE(0.0, tend, 100);
end

