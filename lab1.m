%
% AUTHOR: Zheng Wang 
% UCD ID: 19206109
%
% FUNCTION:
% This program simulate the 8-QAM digital communnication system;
% It plots the simulated and theoretical curve of SER vs Es/E0;
% It also plots the simulated curve of BER vs Es/E0
%
%%
clc;
clear;
close all;

N = 3e5; % Number of bits which is utlized in this assignment
M = 8; % M=8
BPS = log2(M); % Bits per symbol
threshold = 0.5; % set the threshold distinguishing between 0 and 1.

% Bits-to-Symbol mapping using Gray-Coded, which is in the 
% form of compex number accorind to the constellation diagram.
map=[-1+3j, -1+1j, -1-3j, -1-1j, 1+1j, 1+3j, 1-3j, 1-1j];

% plot the constellation diagram for the 8-QAM system.
figure(1)
plot(map, 'r.', 'MarkerSize', 30);
grid on;
title('Constellation Diagram for 8-QAM system');
axis([-4 4 -4 4]);
xlabel('Inphase');
ylabel('Quadrature');

set(gcf,'unit','centimeters','position',[7 5 15 10])
set(gca,'Position',[.15 .15 .75 .75]);
set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',0.4);
set(gca,'FontSize',12,'FontName','monospace','XColor','k','YColor','k','linewidth',2);
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',0.4,'color','w');

%%
EsN0_dB = 0:20; % define the range of EsN0 in dB
EbN0_dB = EsN0_dB - 10*log10(BPS); % define the range of EbN0 in dB

SimulatedBER = zeros(1,length(EbN0_dB));     % simulational values of BER
SimulatedSER = zeros(1,length(EbN0_dB));     % simulational values of SER
theoreticalSER = zeros(1,length(EbN0_dB));   % theoretical value of SER
n = 1;  % counter

% normalize the symbols with its average symbol energy
norm_map = sqrt(1/6) * map;
% get the real components of the normalized mapping
In = real(norm_map);
% get the complex components of the normalized mapping
Qu = imag(norm_map);

for i = EbN0_dB
    % generate input data
    bits = double(rand(1,N)>=threshold); % random bits generation
    input_bits = reshape(bits, BPS, [])'; % reshape the bits 
    input_bits = input_bits*(2.^((BPS-1):-1:0)).'; % convert bits to symble
    symbol = norm_map(input_bits + 1).'; % use index values to map Constellation 
    EbN0 = 10.^(i/10); % generate EbN0 in linear
    EsN0 = 10.^(EsN0_dB(n)/10); % convert EsN0 from dB to linear
    Noise = sqrt(1./(2*BPS*EbN0)); 
    w = Noise*(randn(1,length(symbol))+1i*randn(1,length(symbol)))'; % generate the noises
    y = symbol + w; % generate the final signal
    
    % demodulization Process
    DemoduleS = zeros(1, length(y));
    for j = 1:length(y)
        % obtain the minimum Eucledian Distance to find the demodulized
        % value
        [mValue, mIndex] = min(sqrt((real(y(j))-In(1,:)).^2+(imag(y(j))-Qu(1,:)).^2));
        DemoduleS(j) = mIndex - 1;
    end
    
    % 1. BER
    % convert symbols to bits
    DemoduleBits = dec2bin(DemoduleS)-'0';
    output_bits = reshape(DemoduleBits.', 1, [])'; % reshape the demoulated bits
   
    % Compute simulation BER
    bitErrors = sum(sum(xor(output_bits.', bits))); % find the number of the error bits
    SimulatedBER(n) = bitErrors / N; % calculate simulated BER
    
    % 2. SER
    % 2.1 Compute the theoretical SER according to the derived expression 
    theoreticalSER(n) = (5/2)*qfunc(sqrt(EsN0/3))- ...
    (3/2)*qfunc(sqrt(EsN0/3))*qfunc(sqrt(EsN0/3)); 

    % 2.2 Compute Simulated SER 
    SimulatedSER_temp = input_bits.' - DemoduleS;
    SimulatedSER(n) = sum(SimulatedSER_temp(:)~=0) / (N/BPS);
    
    % Increase index by 1 
    n = n + 1;
    
end

%%
% Plot simulated and theoretical SER curves in the same figure
figure(2);
semilogy(EsN0_dB,SimulatedSER, '-*', 'LineWidth',1.5);
hold on;
semilogy(EsN0_dB,theoreticalSER,'--o', 'LineWidth',1.5);
grid on;
legend('simulation','theoretical'); 
xlabel('E_s/N_0 dB');
ylabel('SER');
title('The symbol error rate (SER) versus Es/N0 and the theoretical SER curve for the system');

set(gcf,'unit','centimeters','position',[7 5 15 10])
set(gca,'Position',[.15 .15 .75 .75]);
set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',0.4);
set(gca,'FontSize',12,'FontName','monospace','XColor','k','YColor','k','linewidth',2);
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',0.4,'color','w');

% Plot simulated BER curve
figure(3);
semilogy(EbN0_dB, SimulatedBER, 'b-*', 'LineWidth',2);
grid on;
title('The bit error rate (BER) versus Eb/No for 8-ary Modulation');
legend('Simulated BER');
xlabel('E_b/N_0 dB');
ylabel('BER');

set(gcf,'unit','centimeters','position',[7 5 15 10])
set(gca,'Position',[.15 .15 .75 .75]);
set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',0.4);
set(gca,'FontSize',12,'FontName','monospace','XColor','k','YColor','k','linewidth',2);
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',0.4,'color','w');







        


