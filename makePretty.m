function makePretty()
% prettify plots

alw = 1.2;    % AxesLineWidth
fsz = 16;      % Fontsize
set(gca, 'FontSize', fsz, 'LineWidth', alw);
%set(gcf, 'FontSize', fsz);
grid on;
end

