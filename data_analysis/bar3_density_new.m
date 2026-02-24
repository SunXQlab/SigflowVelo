close all
clear
clc

data = xlsread('../umap_coords/Dutcal.xlsx', 'Sheet1');
x = data(:,1);
y = data(:,2);
xlist = linspace(min(x), max(x));
ylist = linspace(min(y), max(y));
[c, xmesh1, ymesh1, zmesh1] = density2C(x, y, xlist, ylist, 15, 15);

data2 = xlsread('../umap_coords/Dutcal.xlsx', 'Sheet2');
x2 = data2(:,1);
y2 = data2(:,2);
xlist2 = linspace(min(x2), max(x2));
ylist2 = linspace(min(y2), max(y2));
[c2, xmesh2, ymesh2, zmesh2] = density2C(x2, y2, xlist2, ylist2, 15, 15);

custom_colors = [
    112, 201, 235;   
    202, 232, 242;
    243, 213, 216;
    255, 152, 150;
    233, 093, 105;
];
custom_colors = custom_colors / 255;

num_colors = 100;
custom_colors = interp1(linspace(0, 1, size(custom_colors, 1)), custom_colors, linspace(0, 1, num_colors));

axisFontSize = 14;
titleFontSize = 16;
tickFontSize = 12;
colorbarFontSize = 12;

figure;

subplot(2, 1, 1);
b1 = bar3(zmesh1);
xlabel('UMAP 1', 'FontSize', axisFontSize);
ylabel('UMAP 2', 'FontSize', axisFontSize);
zlabel('Density', 'FontSize', axisFontSize);
colormap(custom_colors);
caxis([0, 0.04]);
zlim([0 0.04]);
for k = 1:length(b1)
    zdata = b1(k).ZData;
    b1(k).CData = zdata;
    b1(k).FaceColor = 'interp';
end
xticks([]);
yticks([]);
view(49, 21);

set(gca, 'FontSize', tickFontSize);

axes_pos = get(gca, 'Position');
text('String', 'Dutcal cells in sensitive group', 'Position', [axes_pos(1) - 0.1, 1.05, 0], 'Units', 'normalized', 'HorizontalAlignment', 'left', 'FontSize', titleFontSize);

subplot(2, 1, 2);
b2 = bar3(zmesh2);
xlabel('UMAP 1', 'FontSize', axisFontSize);
ylabel('UMAP 2', 'FontSize', axisFontSize);
zlabel('Density', 'FontSize', axisFontSize);
colormap(custom_colors);
caxis([0, 0.04]);
zlim([0 0.04]);
for k = 1:length(b2)
    zdata = b2(k).ZData;
    b2(k).CData = zdata;
    b2(k).FaceColor = 'interp';
end
xticks([]);
yticks([]);
view(49, 21);

set(gca, 'FontSize', tickFontSize);

axes_pos2 = get(gca, 'Position');
text('String', 'Dutcal cells in resistant group', 'Position', [axes_pos2(1) - 0.1, 1.05, 0], 'Units', 'normalized', 'HorizontalAlignment', 'left', 'FontSize', titleFontSize);

set(gcf, 'Position', [100, 100, 800, 1000]);

cb = colorbar;
cb.Position = [0.8, 0.25, 0.02, 0.4];
cb.FontSize = colorbarFontSize;

pos1 = get(subplot(2,1,1), 'Position');
pos2 = get(subplot(2,1,2), 'Position');
subplots_adjust = 0.85;
pos1(3) = subplots_adjust * pos1(3);
pos2(3) = subplots_adjust * pos2(3);
set(subplot(2,1,1), 'Position', pos1);
set(subplot(2,1,2), 'Position', pos2);
