%
% AUTHOR: Zheng Wang 
% UCD ID: 19206109
%
% FUNCTION:
% This program simulate the 16-QAM digital communnication system 
% with Hamming code scheme;
% It plots the simulated and theoretical curve of SER vs Es/E0;
% It also plots the simulated curve of BER vs Es/E0;
% Moreover, it compares the simulated BER results between the scenarios
% that when Hamming code is used and when it is not used.

%%
clc;
clear;
close all;

N = 3e5;                    % Number of bits to be generated
M = 16;                     % 16-ary QAM modulation 
BPS = log2(M);              % Bits per symbol

EbN0_dB = -5:20;             % define EbN0 in dB
EbN0 = 10.^(EbN0_dB/10);    % convert EbN0_dB to linear
EsN0 = log2(M)*EbN0;        % define EsN0

SymbolErrors = zeros(1,length(EbN0_dB));  % symbol errors
BER_simulated = zeros(1,length(EbN0_dB));  % simulational value of BER
BER_theoretical = zeros(1,length(EbN0_dB)); % theoretical value of BER
SER_theoretical = zeros(1,length(EbN0_dB)); % theoretical value of SER
SER_simulated = zeros(1,length(EbN0_dB)); % simulational value of SER
BER_coded = zeros(1,length(EbN0_dB)); % simulational value of BER with hamming code
BER_uncoded = zeros(1,length(EbN0_dB)); % theoretical value of BER with hamming code
input_bits = zeros(1, N/BPS);
det = zeros(1, N/log2(M));
dis = zeros(1,M);
%% Detemine the hamming-coded bits
hamming_code = zeros(1, N/4*7);   % 7*4 hamming code, note that N/4*7 is  
                                  % the leng of the data to be encoded
bits = round(rand(1,N)); % The input bits

% Obtain the hamming coded redults  
for i = 0:N/4-1 
    for j = 1:7
        if(j==5)
            hamming_code(1, i*7+5) = xor(bits(1, i*4+2),...
            xor(bits(1, i*4+3), bits(1, i*4+4)));
        elseif(j==6)
            hamming_code(1, i*7+6) = xor(bits(1, i*4+1),...
            xor(bits(1, i*4+2), bits(1, i*4+4)));
        elseif(j==7)
            hamming_code(1, i*7+7) = xor(bits(1, i*4+1),...
            xor(bits(1, i*4+3), bits(1, i*4+4)));
        else
             hamming_code(1, i*7+j) = bits(1, i*4+j);
        end
    end           
end

%% Modulation Process
% Define constellation 
map = [-3+3j, -1+3j, 1+3j, 3+3j, -3+1j, -1+1j, 1+1j, 3+1j,... 
       -3-1j, -1-1j, 1-1j, 3-1j, -3-3j, -1-3j, 1-3j, 3-3j];

% Define Grey-Coded Symbols    
matching = [[0,0,0,0];[0,0,0,1];[1,0,0,1];[1,0,0,0];...
            [0,0,1,0];[0,0,1,1];[1,0,1,1];[1,0,1,0];...
            [0,1,1,0];[0,1,1,1];[1,1,1,1];[1,1,1,0];...
            [0,1,0,0];[0,1,0,1];[1,1,0,1];[1,1,0,0]];
        
% Define Hamming-Code decoder Syndrome
matching2 = [[1,1,1];[1,0,1];[1,1,0];[0,1,1];[1,0,0];[0,1,0];[0,0,1]];

% Mapping
for j = 1: N/BPS
    temp_bits = bits(4*j-3:4*j);
    for k = 1:16
        if(temp_bits == matching(k,:))
            input_bits(j) = map(k);
        end
    end
end

%% 16-QAM demodulation process
N0 = (sum(abs(map).^2) / length(map))./EsN0;
for n = 1:length(EbN0_dB)
    Noise = sqrt(N0(n)/2); % Define the input noise
    w = Noise * randn(1, N/log2(M))+1j*sqrt(N0(n)/2)*randn(1,N/log2(M));
    y = input_bits + w; % generate the final signal 
    for k = 1:N/log2(M)
        for h = 1:M
            dis(h) = norm(y(k) - map(h))^2;
        end
        % obtain the minimum Eucledian Distance to find the demodulized value
        i = find(dis == min(dis));
        det(k) = map(i);
        if(det(k)~=input_bits(k))
            SymbolErrors(n) = SymbolErrors(n)+1;
        end
    end
    
    % Compute the simulated value of SER
    SER_simulated(n) = SymbolErrors(n) / (N/log2(M)); 
    
    % Compute the theoretical value of SER
    Coding_Gain = (2^4-1)*3*2^2/(6*18);  % compute coding gain
    P_bit = ((13/14).*(erfc((sqrt(Coding_Gain)*2)./(2.*sqrt((7/4)./EbN0)))))./2;
    SER_theoretical(n) = 3*qfunc(sqrt(4/5 * EbN0(n))) * (1-3/4*qfunc(sqrt(4/5*EbN0(n))));
    
    % Compute simlulated and theoretical values of BER 
    BER_simulated(n) = (1/log2(M)) * SER_simulated(n);
    BER_theoretical(n) = (1/log2(M)) * SER_theoretical(n);
end

%% Demodulation when hamming code is used 
input_bits0 = (hamming_code - 0.5) * 2;
N1 = (7/4)./EbN0;
for i = 1:length(EbN0_dB)
    Noise = sqrt(N1(i)/2);
    w = Noise * randn(1,N/4*7);
    y = input_bits0 + w;
    Demodule = zeros(1, N/4*7);    % demoudulator container
    Decode = zeros(1, N);          % decoder container
    
    % Demodulate the input signal with (7,4) hamming-coded
    for m = 1 : (N/4*7)
        if(y(m)>=0)
            Demodule(m) = 1;
        else
            Demodule(m) = 0;
        end
    end
    
    % decoding the hamming-coded signal 
    for a = 0 : ((N/4*7)/7-1)
        % The decoded bit-stream
        Decode(1,a*4+1) = Demodule(1,a*7+1);
        Decode(1,a*4+2) = Demodule(1,a*7+2);
        Decode(1,a*4+3) = Demodule(1,a*7+3);
        Decode(1,a*4+4) = Demodule(1,a*7+4);
        
        % The parity check matrix verification
        X(1) = xor(Demodule(1,a*7+1),xor(Demodule(1,a*7+2),...
                xor(Demodule(1,a*7+3),Demodule(1,a*7+5))));
        X(2) = xor(Demodule(1,a*7+1),xor(Demodule(1,a*7+3),...
                xor(Demodule(1,a*7+4),Demodule(1,a*7+6))));
        X(3) = xor(Demodule(1,a*7+1),xor(Demodule(1,a*7+2),...
                xor(Demodule(1,a*7+4),Demodule(1,a*7+7))));
        for b = 1:7
            if(X == matching2(b,:))
                Demodule(1,a*4+b) = ~Demodule(1,a*4+b);
            end
        end
    end
    
    % count the error bit
    for m = 1 : N
        if(Decode(m) ~= bits(m))
            SymbolErrors(i) = SymbolErrors(i) + 1;
        end
    end
    
    % Compute simulated BER
    BER_coded(i) = (SymbolErrors(i)/ N);
    
    % Compute simulated BER without hamming-coded
    BER_undecoded = 1-(1-P_bit).*(1-P_bit);
end

%%
% Plot simulated SER and theoretical SER
figure(1)
semilogy(EbN0_dB, SER_simulated, 'b-*', 'LineWidth', 2);
hold on;
semilogy(EbN0_dB, SER_theoretical, '--o', 'LineWidth', 2);
legend('simulation','theoretical'); 
xlabel('E_s/N_0 dB');
ylabel('SER');
grid on;
axis([0 16 10^-5 1]);
title('The symbol error rate (SER) versus Es/N0 and the theoretical SER curve for the system');

set(gcf,'unit','centimeters','position',[7 5 15 10])
set(gca,'Position',[.15 .15 .75 .75]);
set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',0.4);
set(gca,'FontSize',12,'FontName','monospace','XColor','k','YColor','k','linewidth',2);
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',0.4,'color','w');

% Plot simulated BER and theoretical BER
figure(2)
Fig4 = semilogy(EbN0_dB,BER_simulated,'b-*','LineWidth',2);
hold on
Fig5 = semilogy(EbN0_dB,BER_theoretical,'--o','LineWidth',2);
legend([Fig4(1),Fig5(1)],'16-QAM simulation BER','16-QAM theoretical BER');
xlabel('E_b/N_0 dB');
ylabel('BER');
grid on
axis([0 14 10^-6 1]);
title('The bit error rate (BER) versus Eb/N0 curve for the system');

set(gcf,'unit','centimeters','position',[7 5 15 10])
set(gca,'Position',[.15 .15 .75 .75]);
set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',0.4);
set(gca,'FontSize',12,'FontName','monospace','XColor','k','YColor','k','linewidth',2);
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',0.4,'color','w');

% Plot the comparison of the simulated BER results between the scenarios
% that when Hamming code is used and when it is not used.
figure(3)
Fig1 = semilogy(EbN0_dB,BER_coded,'r-','LineWidth',3.5);
hold on
Fig2 = semilogy(EbN0_dB,BER_undecoded,'b-','LineWidth',3.5);
legend([Fig1(1),Fig2(1)],'16-QAM with (7,4)hamming code','16-QAM without (7,4)hamming code'); 
axis([-6 20 10^-5 1]); 
xlabel('E_s/N_0 dB');
ylabel('BER');
grid on
title('The BER versus Es/N0 with hamming code and that without hammingcode curves for the system');

set(gcf,'unit','centimeters','position',[7 5 15 10])
set(gca,'Position',[.15 .15 .75 .75]);
set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',0.4);
set(gca,'FontSize',12,'FontName','monospace','XColor','k','YColor','k','linewidth',2);
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',0.4,'color','w');
        
        
            
             
                
        
    
    

    
            
        
    
    
    

    
        
    
    



