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


figure(1)
plot(restart, (gmres_gs_time./gmres_time)*100, lineSpec{1}, 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor',  markerFill{1});
grid on
hold on
plot(restart, (timeGSAll(1,:)./timeAll(1, :))*100, lineSpec{2}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor',  markerFill{2});
plot(restart, (timeGSAll(2,:)./timeAll(2, :))*100, lineSpec{3}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor',  markerFill{3});
plot(restart, (timeGSAll(3,:)./timeAll(3, :))*100, lineSpec{4}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor',  markerFill{4});
plot(restart, (timeGSAll(4,:)./timeAll(4, :))*100, lineSpec{5}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor',  markerFill{5});
ax=gca;
ax.FontSize=18;

xlabel('restart','FontSize', 18, 'FontWeight','bold')
ylabel('ratio: time for GS to time to solution (%)', 'FontSize', 18, 'FontWeight','bold')
ll1 = legend( 'GMRES-MGS (HYPRE)', 'GMRES-MGS (new implementation)', 'GMRES-CGS1 w/alt norm',...
    'GMRES-CGS2', 'GMRES-two synch', 'Location', 'southeast')
ll1.FontSize=14;
xlim([5 72])
ylim([0 100])


print(gcf, '-djpeg','-r300', 'sc18fig1.jpg');


figure(2)

plot(restart, iters_gmres, lineSpec{1}, 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', markerFill{1});
grid on
hold on
plot(restart, itersAll(1,:), lineSpec{2}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', markerFill{2});
plot(restart, itersAll(2,:), lineSpec{3}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', markerFill{3});
plot(restart, itersAll(3,:), lineSpec{4}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', markerFill{4});
plot(restart, itersAll(4,:), lineSpec{5}, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', markerFill{5});
ax=gca;
ax.FontSize=18;

title('Number of iterations', 'FontSize', 20, 'FontWeight', 'bold');
xlabel ('restart',  'FontSize', 18, 'FontWeight','bold');
ylabel('iterations',  'FontSize', 18, 'FontWeight','bold');

ll2=legend( 'GMRES-MGS (HYPRE)', 'GMRES-MGS (new implementation)', 'GMRES-CGS1 w/alt norm', 'GMRES-CGS2', 'GMRES-two synch')
ll2.FontSize=14;
xlim([5 72])
ylim([0 12000])

print(gcf, '-djpeg','-r300', 'sc18fig3.jpg');
