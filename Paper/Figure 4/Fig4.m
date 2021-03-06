close all;
clear all;

restart = [ 5 10 15 20 25 30 35 40 45 50 55 60 65 70 72 ];
gmres_time  =[20.6186  20.2074  29.8987  21.5269  15.88 31.15 37.7073 60.2905 43.29 71.17  78.1694  87.212 95.4134  93.1558  103.283];
gmres_gs_time = [9.25  13.11  21.64  16.4919 12.87 25.8072 32.33 52.77 38.4927 57.27 63.4808 80.03 88.16 87.56 96.82];
iters_gmres = [ 9346 4233 3690 2283 1421 2088  2143  2366  1656 2004 2012 1970 1917 1756 1749];

load('timeAll.mat');
load('timeGSAll.mat');
load('itersAll');


lineSpec = {'b-o', 'k-<', 'r-d', 'g-s', 'k-d'};
markerFill = {'c', 'g', 'y', 'k', 'r'};
comboData = [ gmres_time./timeAll(1, :); gmres_time./timeAll(2, :); gmres_time./timeAll(3, :);  gmres_time./timeAll(4, :)  ];

figure(1)
plot(restart, gmres_time, lineSpec{1}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', markerFill{1});
grid on
hold on

plot(restart, timeAll(1, :), lineSpec{2}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor',  markerFill{2});
plot(restart, timeAll(2, :), lineSpec{3}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor',  markerFill{3});
plot(restart,timeAll(3, :), lineSpec{4}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor',  markerFill{4});
plot(restart,timeAll(4, :), lineSpec{5}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor',  markerFill{5});
ax=gca;
ax.FontSize=18;

xlabel('restart', 'FontSize', 18, 'FontWeight','bold')
ylabel('time to solution(s)', 'FontSize', 18, 'FontWeight','bold')

ll1 =legend( 'GMRES-MGS (HYPRE)',...
    'GMRES-MGS (new implementation)',...
    'GMRES-CGS1 w/alt norm',...
    'GMRES-CGS2',...
    'GMRES-two synch', 'Location', 'northwest');
ll1.FontSize=14;
xlim([5 72])
ylim([0 110])
print(gcf, '-djpeg','-r300', 'sc18fig0.jpg');

figure(2) 

plot(restart, comboData(1,:), lineSpec{2}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', markerFill{2});
grid on
hold on
plot(restart,  comboData(2,:), lineSpec{3}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', markerFill{3});
plot(restart,  comboData(3,:), lineSpec{4}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', markerFill{4});
plot(restart,  comboData(4,:), lineSpec{5}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', markerFill{5});
  
ax=gca;
ax.FontSize=18;  
grid on
xlabel ('restart', 'FontSize', 18, 'FontWeight','bold');
ylabel('speedup',  'FontSize', 18, 'FontWeight','bold');
title('Speedup in relation to GMRES-MGS (HYPRE)',  'FontSize', 20, 'FontWeight','bold');
ll2=legend(  'GMRES-MGS (new implementation)', 'GMRES-CGS1 w/alt norm', 'GMRES-CGS2', 'GMRES-two synch',...
    'Location', 'northwest');
ll2.FontSize=14;
xlim([5 72])
print(gcf, '-djpeg','-r300', 'sc18fig2.jpg');
