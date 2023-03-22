close
clear all
clc

% https://de.mathworks.com/matlabcentral/answers/102646-how-do-i-add-a-logo-image-into-a-plot-or-a-figure
%% creating fig
screensize = get(0, 'MonitorPositions');
screensize = screensize(1, 3:4);
height = screensize(2)*4.5/5;
width = screensize(1)*4/5;
figPosition = [screensize(2)*1/5 screensize(1)*0.15/5 width height];

% delete old stuff
fig = findall(0, 'tag', 'XAI');
delete(fig);

h.GUI.fig = figure('tag', 'XAI', 'resize', 'off', 'NumberTitle', 'off',...
    'Name', 'XAI', 'Position', figPosition, 'ToolBar', 'None', 'MenuBar', 'None');

try
    frame_h = get(h.GUI.fig, 'javaframe');
    frame_h.setFigureIcon(javax.swing.ImageIcon('Images/UniUlm.png'));
catch
end

%% Load Data
try
    table = readtable('input/train.csv');
    i = 1;
    while true
        h.data{1, i} = readtable(['explanations/explanations_Store' num2str(i) '.csv']);
        str{i} = ['Store ' num2str(i)];
        index = find(table.Store == i & table.Sales ~= 0);
        h.data{2,i} = table(index, [3,4]);
        i = i+1;
    end
catch 
end
%% creating GUI
% title
uicontrol(h.GUI.fig, 'style', 'text', 'String', 'Rossmann Store Sales Forecasting with XAI',...
    'units', 'normalized', 'Position', [0.11 0.9 0.5, 0.05], 'HorizontalAlignment', 'left',...
    'FontSize', 20, 'FontWeight', 'bold');

% dropdownmenu
uicontrol(h.GUI.fig, 'style', 'Popupmenu', 'String', str,...
    'units', 'normalized', 'Position', [0.125 0.87 0.1, 0.025], 'HorizontalAlignment', 'left',...
    'FontSize', 10, 'Callback', @callbackChangeStore);

% plot
h.GUI.Plot = axes(h.GUI.fig, 'position', [0.05 0.37 0.7 0.45]);
% predicted data
X_DataP = table2array(h.data{1,1}(:, 6));
Y_DataP = table2array(h.data{1,1}(:, 5));
Z_DataP = h.data{1, 1}.Anchors;
index = find(h.data{2,1}.Date == X_DataP(end));
% train dat
X_DataT = table2array(h.data{2,1}(index:end, 1));
Y_DataT = table2array(h.data{2,1}(index:end, 2));
% plot data
plot(h.GUI.Plot, X_DataT, Y_DataT, 'b');
hold on
plot(h.GUI.Plot, X_DataP, Y_DataP, 'r');
h.GUI.Plot.XLim   = [X_DataP(1)-calmonths(8) X_DataP(1)];

% Ereignis-Handler für jeden Punkt hinzufügen
for i = 1:length(X_DataP)
    h.GUI.Scatter(i) = scatter(h.GUI.Plot, X_DataP(i), Y_DataP(i), 5, 'filled', 'r');
    set(h.GUI.Scatter(i), 'ButtonDownFcn', {@pointClickCallback, Z_DataP{i}, h.data{1, 1}.Precisions(i), h.data{1, 1}.PredictedSales(i)});
end
legend('Trained Sales Data', 'Predicted Sales Data');
% By default the first predicted point is selected
h.GUI.Scatter(end).CData = [0 1 0];
h.GUI.Scatter(end).SizeData = 40;
h.Marked = h.GUI.Scatter(end).SeriesIndex;
% some general guidelines
FS_H = 13;
FS_SH = 11;
FS_N = 10;
height = 0.025;

%% Explanation for default point
% exact date
h.GUI.Date = uicontrol(h.GUI.fig, 'style', 'text', 'String', datestr(X_DataP(end)),...
    'units', 'normalized', 'Position', [0.77 0.77 0.5, 0.05], 'HorizontalAlignment', 'left',...
    'FontSize', FS_H, 'FontWeight', 'bold');
% number of predicted Sales
uicontrol(h.GUI.fig, 'style', 'text', 'String', 'Number of Predicted Sales by AI:',...
    'units', 'normalized', 'Position', [0.78  0.725 0.5, 0.05], 'HorizontalAlignment', 'left',...
    'FontSize', FS_SH, 'FontWeight', 'bold');
h.GUI.Pred = uicontrol(h.GUI.fig, 'style', 'text', 'String', round(h.data{1,1}.PredictedSales(end)),...
    'units', 'normalized', 'Position', [0.785 0.725-height 0.5, 0.05], 'HorizontalAlignment', 'left',...
    'FontSize', FS_N);
% title explanation
uicontrol(h.GUI.fig, 'style', 'text', 'String', 'Explanation:',...
    'units', 'normalized', 'Position', [0.78 0.65 0.5, 0.05], 'HorizontalAlignment', 'left',...
    'FontSize', FS_SH, 'FontWeight', 'bold');
% split the explantation
strs = strsplit(h.data{1, 1}.Anchors{end}, ',');
% plot every point of explanation
for i=1:length(strs)
    h.GUI.Exp{i} = uicontrol(h.GUI.fig, 'style', 'text', 'String', ['- ' strs{i}(3:end-1)],...
    'units', 'normalized', 'Position', [0.78 0.65-i*height 0.5, 0.05], 'HorizontalAlignment', 'left',...
    'FontSize', FS_N);
end

% Precision of Explanation
h.GUI.TitlePrec = uicontrol(h.GUI.fig, 'style', 'text', 'String', 'Precision:',...
    'units', 'normalized', 'Position', [0.78  0.65-(i+2)*height 0.5, 0.05], 'HorizontalAlignment', 'left',...
    'FontSize', FS_SH, 'FontWeight', 'bold');
h.GUI.Prec = uicontrol(h.GUI.fig, 'style', 'text', 'String', h.data{1,1}.Precisions(end),...
    'units', 'normalized', 'Position', [0.785 0.65-(i+3)*height 0.5, 0.05], 'HorizontalAlignment', 'left',...
    'FontSize', FS_N);

%% General Informations about the Store
% Specific Store and All Stores
uicontrol(h.GUI.fig, 'style', 'text', 'String', 'Mean All Stores:',...
    'units', 'normalized', 'Position', [0.35 0.27 0.5, 0.05], 'HorizontalAlignment', 'left',...
    'FontSize', FS_H, 'FontWeight', 'bold');

h.GUI.GI.NR = uicontrol(h.GUI.fig, 'style', 'text', 'String', 'Store 1:',...
    'units', 'normalized', 'Position', [0.5 0.27 0.5, 0.05], 'HorizontalAlignment', 'left',...
    'FontSize', FS_H, 'FontWeight', 'bold');

str = {'Average Sales Train Data:', 'Average Sales Predicted Data:', 'Standard Deviation of Predicted Sales:',...
    'Average Precision of Explanation:'};
dataS = {mean(Y_DataT), mean(Y_DataP), sqrt(sum((h.data{1, 1}.ESales-Y_DataP).^2)/length(Y_DataP)), mean(h.data{1,1}.Precisions)};

for i = 1:size(h.data, 2)
    meanST(i) = mean(h.data{2, i}.Sales);
    meanSP(i) = mean(h.data{1, i}.PredictedSales);
    Dev(i) = sqrt(sum((h.data{1, i}.ESales-h.data{1, i}.PredictedSales).^2)/length(h.data{1, i}.PredictedSales));
    meanPrec(i) = mean(h.data{1, i}.Precisions);
end
dataAS = {mean(meanST), mean(meanSP), mean(Dev), mean(meanPrec)};

height = 0.03;
for i = 1:length(str)
    uicontrol(h.GUI.fig, 'style', 'text', 'String', str{i},...
        'units', 'normalized', 'Position', [0.125 0.225-height*(i-1) 0.2, 0.05], 'HorizontalAlignment', 'left',...
        'FontSize', FS_SH, 'FontWeight', 'bold');
    uicontrol(h.GUI.fig, 'style', 'text', 'String', dataAS{i},...
        'units', 'normalized', 'Position', [0.35 0.22-height*(i-1) 0.5, 0.05], 'HorizontalAlignment', 'left',...
        'FontSize', FS_N);
    h.GUI.GI.Data{i} = uicontrol(h.GUI.fig, 'style', 'text', 'String', dataS{i},...
        'units', 'normalized', 'Position', [0.5 0.22-height*(i-1) 0.5, 0.05], 'HorizontalAlignment', 'left',...
        'FontSize', FS_N);
end

% Current Store
h.StoreIdx = 1;
% store data
guidata(h.GUI.fig, h);
%% Funktion for XAI
function pointClickCallback(src, event, strs, prec, pred)
    % get data
    h = guidata(src);
    % undo old Marker
    h.GUI.Scatter(h.Marked-2).CData = [1 0 0];
    h.GUI.Scatter(h.Marked-2).SizeData = 5;
    % get the recent point
    x = get(src, 'XData');
    y = get(src, 'YData');
    % new Marker
    src.SizeData = 40;
    src.CData = [0 1 0];
    h.Marked = src.SeriesIndex;
    
    % new Prediction
    h.GUI.Pred.String = round(pred);
    
    % Explanation
    strs = strsplit(strs, ',');
    h.GUI.Date.String = datestr(x);

    height = 0.025;
    for i=1:length(strs)
        h.GUI.Exp{i}.String = ['- ' strs{i}(3:end-1)];
    end
    for j = length(strs)+1:length(h.GUI.Exp)
        h.GUI.Exp{j}.String = '';
    end
    
    h.GUI.TitlePrec.Position = [0.78  0.65-(i+2)*height 0.5, 0.05];
    h.GUI.Prec.String = prec;
    h.GUI.Prec.Position = [0.785  0.65-(i+3)*height 0.5, 0.05];
    
    % store data
    guidata(src, h);
end

%% function for changing store
function callbackChangeStore(src,event)
     % get data
    h = guidata(src);
    % new data
    %% General Informations about the Store
    % Specific Store and All Stores
    idx = src.Value;
    h.GUI.GI.NR.String =  [src.String{idx} ':'];

  
    dataS = {mean(h.data{2, idx}.Sales), mean(h.data{1, idx}.PredictedSales),...
        sqrt(sum((h.data{1, idx}.ESales-h.data{1, idx}.PredictedSales).^2)/length(h.data{1, idx}.PredictedSales)), ...
        mean(h.data{1,idx}.Precisions)};


    for i = 1:length(dataS)
        h.GUI.GI.Data{i}.String = dataS{i};
    end
    
    %% plot new data
    cla(h.GUI.Plot);
    % predicted data
    X_DataP = table2array(h.data{1,idx}(:, 6));
    Y_DataP = table2array(h.data{1,idx}(:, 5));
    Z_DataP = h.data{1, idx}.Anchors;
    index = find(h.data{2,idx}.Date == X_DataP(end));
    % train dat
    X_DataT = table2array(h.data{2,idx}(index:end, 1));
    Y_DataT = table2array(h.data{2,idx}(index:end, 2));
    % plot data
    plot(h.GUI.Plot, X_DataT, Y_DataT, 'b');
    hold on
    plot(h.GUI.Plot, X_DataP, Y_DataP, 'r');
    h.GUI.Plot.XLim   = [X_DataP(1)-calmonths(8) X_DataP(1)];

    % Ereignis-Handler für jeden Punkt hinzufügen
    for i = 1:length(X_DataP)
        h.GUI.Scatter(i) = scatter(h.GUI.Plot, X_DataP(i), Y_DataP(i), 5, 'filled', 'r');
        set(h.GUI.Scatter(i), 'ButtonDownFcn', {@pointClickCallback, Z_DataP{i}, h.data{1, idx}.Precisions(i), h.data{1, idx}.PredictedSales(i)});
    end
    legend('Trained Sales Data', 'Predicted Sales Data');
    % By default the first predicted point is selected
    try
        h.GUI.Scatter(h.Marked-2).CData = [1 0 0];
        h.GUI.Scatter(h.Marked-2).SizeData = 5;
    catch
    end
    
    h.GUI.Scatter(i).CData = [0 1 0];
    h.GUI.Scatter(i).SizeData = 40;
    h.Marked = h.GUI.Scatter(i).SeriesIndex;
    
    %% Explanation
    %% Explanation for default point
    % exact date
    h.GUI.Date.String = datestr(X_DataP(end));
    % number of predicted Sales
    h.GUI.Pred.String = round(h.data{1,idx}.PredictedSales(end));
    
    
    strs = strsplit(h.data{1, idx}.Anchors{end}, ',');
    % plot every point of explanation
    height = 0.025;
    for i=1:length(strs)
        h.GUI.Exp{i}.String = ['- ' strs{i}(3:end-1)];
    end
    for j = length(strs)+1:length(h.GUI.Exp)
        h.GUI.Exp{j}.String = '';
    end
    
    h.GUI.TitlePrec.Position = [0.78  0.65-(i+2)*height 0.5, 0.05];
    h.GUI.Prec.String = h.data{1,idx}.Precisions(end);
    h.GUI.Prec.Position = [0.785  0.65-(i+3)*height 0.5, 0.05];
    
    % Current Store
    h.StoreIdx = idx;
    % store data
    guidata(src, h);
end