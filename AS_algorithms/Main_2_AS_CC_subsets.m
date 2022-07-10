%% Main_2_AS_CC_subsets = Channel capacity VS number of selected subset


%% REF PAPER
% [1] X. Zhou, B. Bai and W. Chen, "Invited Paper: Antenna selection in energy efficient MIMO systems: A survey"
% [2] Gharavi-Alkhansari M, Greshman A. "Fast antenna selection in MIMO systems."

% ----------------------------- MATLAB System -----------------------------
clc
clear;
close all;
format short;
rng('shuffle');
warning('off');

% Antennas parameters
Nr=16;
Nt=16;

% ------------------- SNR parameters -------------------
SNRdB=20;
SNR= 10^(SNRdB/10);

simulation=300;

% Channel capacity butter
capacityOfExAver=[];
capacityOfNBSAver=[];
capacityOfFastAver=[];
capacityOfRandomAver=[];

% Lr = Nr/2
for Lr=1:(10)
    fprintf('Nt(%d), Nr(%d), Lr(%d), SNR(%d)... \n',Nt,Nr,Lr,SNRdB);
    
    % Optimal-subset selection (Conventional scheme)
    capacityOfExSum=0;
    capacityOfNBSSum=0;
    capacityOfFastSum=0;
    capacityOfRandomSum=0;
    
    for sim=1:simulation
        
        % Rayleigh Channel (Channel matrix H)
        % Normalization of channel matrix (dataset version)
        H = abs(sqrt(1/2)*(randn(Nr,Nt)+1j*randn(Nr,Nt)));
        hMin = min(H,[],2);
        hMax = max(H,[],2);
        H = (H-hMin)./(hMax-hMin);
        
        % Number of maximum antenna
        fullAntenna=(1:Nr);
        
        % for Exhaustive search
        TSSB = tic;
        antennaSubset=nchoosek((1:Nr),Lr);
        TFSB = toc(TSSB);
        
        TSF_EX = tic;
        % Exhaustive search
        [capacityOfExSelected]=AS_Exhaustive(Nr,Nt,Lr,SNR,H,antennaSubset);
        capacityOfExSum = capacityOfExSum+capacityOfExSelected;
        TFF_EX = toc(TSF_EX);
        
        TSF_NBS = tic;
        % NBS, Ref [1]
        [capacityOfNBSSelected]=AS_NBS(Nr,Nt,Lr,SNR,H,fullAntenna);
        capacityOfNBSSum=capacityOfNBSSum+capacityOfNBSSelected;
        TFF_NBS = toc(TSF_NBS);
        
        TSF_PF = tic;
        % Fast, Ref [2]
        [capacityOfFastSelected]=AS_Fast(Nr,Nt,Lr,SNR,H,fullAntenna);
        capacityOfFastSum=capacityOfFastSum+capacityOfFastSelected;
        TFF_PF = toc(TSF_PF);
        
        TSF_R = tic;
        % Random
        [capacityOfRandomSelected]=AS_Ran(Nr,Nt,Lr,SNR,H,fullAntenna);
        capacityOfRandomSum=capacityOfRandomSum+capacityOfRandomSelected;
        TFF_R = toc(TSF_R);
        
    end
    TSC_EX = tic;
    capacityOfExAver = [capacityOfExAver, capacityOfExSum/simulation];
    TFC_EX = toc(TSC_EX);
    
    TSC_NBS = tic;
    capacityOfNBSAver=[capacityOfNBSAver,capacityOfNBSSum/simulation];
    TFC_NBS = toc(TSC_NBS);
    
    TSC_PF = tic;
    capacityOfFastAver=[capacityOfFastAver,capacityOfFastSum/simulation];
    TFC_PF = toc(TSC_PF);
    
    TSC_R = tic;
    capacityOfRandomAver=[capacityOfRandomAver,capacityOfRandomSum/simulation];
    TFC_R = toc(TSC_R);
end
% calculation of each antenna selection schemes
% 10^3 -> Decimal removal
T_EX = (TFSB + TFF_EX + TFC_EX)*10^3;
T_NBS = (TFF_NBS + TFC_NBS)*10^3;
T_PF = (TFF_PF + TFC_PF)*10^3;
T_R = (TFF_R + TFC_R)*10^3;

% figure font
figure1 = figure('Color','white');
axes1 = axes('Parent',figure1,'FontName','Times New Roman');
box(axes1,'on'); grid(axes1,'on'); hold(axes1,'all');

Xaxis=(1 : Lr);
plot (Xaxis, capacityOfExAver,'-*','LineWidth', 1.5, 'MarkerSize', 7.5, 'color', '#469B4E'); hold on
plot(Xaxis, capacityOfNBSAver,'-d','LineWidth', 1.5, 'MarkerSize', 7.5, 'color', '#08519C'); hold on 
plot (Xaxis, capacityOfFastAver,'-o','LineWidth', 1.5, 'MarkerSize', 7.5, 'color', '#9ECAE1'); hold on 
plot (Xaxis, capacityOfRandomAver,'-s','LineWidth', 1.5, 'MarkerSize', 7.5, 'color', '#808080'); hold on

% figure set
legend({'Exhaustive','Ref [12]', 'Ref [11]','Random'},'Location','northwest');
title(legend, 'AS schemes', 'FontSize', 11);
set(legend, 'FontName', 'Times New Roman', 'FontSize',11); set(gcf,'Color','w')

% title,x,y labeling
xlabel(['\fontname{times new roman}' 'The number of selected antenna subset'], 'fontsize', 13)
ylabel(['\fontname{times new roman}' 'Channel capacity [bit/s/Hz]'], 'fontsize', 13)

% box set : equal square
axis equal square

% selected subset antenna matrix at least 1
xlim([1 Lr])

% figure box LineWidth
h = gca;
h.LineWidth = 1.15;
hold on; grid on;

%% ------------------------------------------------------------------------------
% figure font
figure2 = figure('Color','white');
axes2 = axes('Parent',figure2,'FontName','Times New Roman');
box(axes2,'on'); grid(axes2,'on'); hold(axes2,'all');

% 0303 legend({'Exhaustive','Ref [12]', 'Ref [11]','Random'},'Location','northwest');
X=categorical({'Exhaustive','Ref [1] ','Ref [2]', 'Random'});
X=reordercats(X,{'Exhaustive','Ref [1]','Ref [2]', 'Random'});
Y=[T_EX, T_NBS,T_PF,T_R];
b=bar(X,Y,'Facecolor','flat', 'FaceAlpha', 0.88, 'FaceColorMode', 'manual', 'LineWidth', 1.05);

xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'FontName','Times New Roman',...
    'fontsize', 9,'HorizontalAlignment','center', 'VerticalAlignment','bottom')

% User-defined hex color code
rgb1  = hex2rgb('#469B4E');
rgb2  = hex2rgb('#08519C');
rgb3  = hex2rgb('#9ECAE1');
rgb4  = hex2rgb('#808080');

% Apply to bar-graph
b.CData(1,:) = [rgb1];    % Greedy Search
b.CData(2,:) = [rgb2];    % Norm-based antenna selection
b.CData(3,:) = [rgb3];    % Fast antenna selection
b.CData(4,:) = [rgb4];    % Random antenna selection

ylim([0 1])

xlabel(['\fontname{times new roman}' 'Antenna selection schemes'], 'fontsize', 12)
ylabel(['\fontname{times new roman}' 'Computation time (sec)'], 'fontsize', 12)

% box set : equal square
axis padded

% figure box LineWidth
h = gca;
h.LineWidth = 1.15;
hold on; grid on;

% add a zoomed zone
% zp = F_PlotZoom();
% zp.plot;

% for draw my favorite color
function [ rgb ] = hex2rgb(hex,range)

assert(nargin>0&nargin<3,'hex2rgb function must have one or two inputs.')

if nargin==2
    assert(isscalar(range)==1,'Range must be a scalar, either "1" to scale from 0 to 1 or "256" to scale from 0 to 255.')
end

if iscell(hex)
    assert(isvector(hex)==1,'Unexpected dimensions of input hex values.')
    if isrow(hex)
        hex = hex';
    end
    hex = cell2mat(hex);
end

if strcmpi(hex(1,1),'#')
    hex(:,1) = [];
end

if nargin == 1
    range = 1;
end
switch range
    case 1
        rgb = reshape(sscanf(hex.','%2x'),3,[]).'/255;
    case {255,256}
        rgb = reshape(sscanf(hex.','%2x'),3,[]).';
    otherwise
        error('Range must be either "1" to scale from 0 to 1 or "256" to scale from 0 to 255.')
end
end