%% Main_1_AS_CC_SNR = Channel capacity VS Signal-to-noise ratio


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
Lr = 14;

SNRMAX = 20;

simulation=5;

% Channel capacity butter
capacityOfExAver=[];
capacityOfNBSAver=[];
capacityOfFastAver=[];
capacityOfRandomAver=[];

for SNRdB=0:SNRMAX
    fprintf('Nt(%d), Nr(%d), Lr(%d), SNR(%d)...\n',Nt,Nr,Lr,SNRdB);
    SNR= 10^(SNRdB/10);
    
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
        antennaSubset=nchoosek((1:Nr),Lr);
        
        % Exhaustive search
        [capacityOfExSelected]=AS_Exhaustive(Nr,Nt,Lr,SNR,H,antennaSubset);
        capacityOfExSum = capacityOfExSum+capacityOfExSelected;

        % NBS, Ref[1]
        [capacityOfNBSSelected]=AS_NBS(Nr,Nt,Lr,SNR,H,fullAntenna);
        capacityOfNBSSum=capacityOfNBSSum+capacityOfNBSSelected;
        
        % Fast, Ref [2]
        [capacityOfFastSelected]=AS_Fast(Nr,Nt,Lr,SNR,H,fullAntenna);
        capacityOfFastSum=capacityOfFastSum+capacityOfFastSelected;
        
        % Random
        [capacityOfRandomSelected]=AS_Ran(Nr,Nt,Lr,SNR,H,fullAntenna);
        capacityOfRandomSum=capacityOfRandomSum+capacityOfRandomSelected;
        
    end
    capacityOfExAver = [capacityOfExAver, capacityOfExSum/simulation];
    capacityOfNBSAver=[capacityOfNBSAver,capacityOfNBSSum/simulation];
    capacityOfFastAver=[capacityOfFastAver,capacityOfFastSum/simulation];
    capacityOfRandomAver=[capacityOfRandomAver,capacityOfRandomSum/simulation];
end

% figure font
figure1 = figure('Color','white');
axes1 = axes('Parent',figure1,'FontName','Times New Roman');
box(axes1,'on'); grid(axes1,'on'); hold(axes1,'all');

Xaxis=(0:SNRMAX);
plot (Xaxis, capacityOfExAver,'-*','LineWidth', 1.5, 'MarkerSize', 7.5, 'color', '#469B4E'); hold on
plot(Xaxis, capacityOfNBSAver,'-d','LineWidth', 1.5, 'MarkerSize', 7.5, 'color', '#08519C'); hold on
plot (Xaxis, capacityOfFastAver,'-o','LineWidth', 1.5, 'MarkerSize', 7.5, 'color', '#9ECAE1'); hold on
plot (Xaxis, capacityOfRandomAver,'-s','LineWidth', 1.5, 'MarkerSize', 7.5, 'color', '#808080'); hold on

% figure set
legend({'Exhaustive','Ref [12]', 'Ref [11]','Random'},'Location','northwest');
title(legend, 'AS schemes', 'FontSize', 11);
set(legend, 'FontName', 'Times New Roman', 'FontSize',11); set(gcf,'Color','w')

% title,x,y labeling
xlabel(['\fontname{times new roman}' 'Signal-to-noise ratio [dB]'], 'fontsize', 13)
ylabel(['\fontname{times new roman}' 'Channel capacity [bit/s/Hz]'], 'fontsize', 13)

% box set : equal square
axis equal square

% figure box LineWidth
h = gca;
h.LineWidth = 1.15;
hold on; grid on;

% add a zoomed zone
%zp = F_PlotZoom();
%zp.plot;