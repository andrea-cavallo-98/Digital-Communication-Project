%
% Andrea Cavallo
% matricola 245715
%
%
% PROJECT #3
% This function plots the eye diagram
%



function [lines] = eyePlot(x,Ntraces,SpT,Offset,Period)
time_plot = linspace(-Period/2,Period/2,SpT);
x = x(Offset:end,:);
figure()
lines = zeros(Ntraces,size((1:SpT),2));
for n = 1:Ntraces
    idx = (n-1)*SpT+(1:SpT);
    lines(n,:) = x(idx);
    plot(time_plot,x(idx), 'b-');
    hold on;
end
end
