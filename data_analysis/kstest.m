clc;
clear;

% 使用 dir 函数列出文件夹下的所有文件和文件夹
files = dir('./umap_coords/');
files(1:2) = [];

celltype_list = {'B cells', 'Dutcal', 'Endothelial', 'Fibroblast', 'Plasma', 'Macrophage', 'Mast', 'NK cells', 'Neutroohils', 'Acinar', 'T cells'};

matrix = [];
% 遍历文件夹中的所有文件和文件夹
for i = 1:length(files)

    filename_list{i} = strcat('./umap_coords/', '', files(i).name);

    sensitive = readmatrix(strcat('./umap_coords/', '', files(i).name), 'Sheet', 1);
    resist = readmatrix(strcat('./umap_coords/', '', files(i).name), 'Sheet', 2);

    x1 = sensitive(:,1);
    y1 = sensitive(:,2);
    x2 = resist(:,1);
    y2 = resist(:,2);

    sen_vector = [x1; y1];
    res_vector = [x2; y2];
    
    % 进行Kolmogorov-Smirnov (K-S) 测试
    [h, p, ks2stat] = kstest2(sen_vector, res_vector);
    
    % 输出结果
    fprintf(strcat(files(i).name), '', 'K-S测试结果:\n');
    fprintf('是否拒绝原假设（两个分布相同）: %d\n', h);
    fprintf('p值: %g\n', p);
    fprintf('统计量值: %g\n', ks2stat);

    
end