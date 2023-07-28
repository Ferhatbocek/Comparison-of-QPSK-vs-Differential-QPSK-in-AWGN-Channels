%% pi/4-DQPSK 

clearvars, clc, close all
N = 20000; % number of bits or symbols
M = 4; % QPSK constellation
k = log2(M); % number of bits per symbol

Eb_N0_dB = 0:15; % multiple Eb/N0 values
SNR_dB = Eb_N0_dB + 10*log10(k); % multiple Es/N0 values

fc = 10.^3; % carrier freq
Tb = 1/fc; % bit duration
t = 0:(Tb/99):Tb;


% Applying pi/4-DQPSK modulation and demodulatin in a single loop 
 
for m = 1:length(Eb_N0_dB)

    % generating random binary sequence	in logical type
    bin_data = rand(1,N)>0.5; % generating 0,1 with equal probability

    % grouping to form of QPSK symbols
    
    Bit_Reshape = reshape(bin_data,2,N/2).'; % Two sequences for binary type
    bintoDecConv = ones(N/2,1)*2.^(k-1:-1:0);	
    Bit_Dec     = sum(Bit_Reshape.*bintoDecConv,2);

    % converting binary to gray coded symbols
    Bit_Gray    = bitxor(Bit_Dec,floor(Bit_Dec/2));
    % generating phase coefficients or multipliers
    Phase_Gray    = 2*Bit_Gray.'+1; 

    % generating differential modulated symbols
    % phi[k] = phi[k-1] + Dphi[k]
    diffPhase = filter(1,[1 -1],Phase_Gray); % start with 0 phase 
    % generating pi/4 DQPSK modulated signal
    dqpsk_signal = exp(1j*diffPhase*pi/4);
    
    % Output of AWGN channel
    y_dqpsk = awgn(dqpsk_signal,Eb_N0_dB(m),'measured');

    % non-coherent demodulation
    % estimated phase information at the receiver
    estPhase  = angle(y_dqpsk);
    % Dphi[k] = phi[k] - phi[k-1]
    est_diffPhase = filter([1 -1],1,estPhase)*4/pi;
    quant_diffPhase = 2*floor(est_diffPhase/2) + 1;  % quantizing
   
    % gray to binary transformation
    quant_diffPhase((quant_diffPhase<0)) = quant_diffPhase((quant_diffPhase<0)) + 8;
    bin_diffPhase     = floor(bitxor(quant_diffPhase,floor(quant_diffPhase/2))/2);
    % estimated binary data after demodulation
    estBit_demod     = (dec2bin(bin_diffPhase.')).'; 
    estBit_demod     = str2num(estBit_demod(1:end).').';

    % counting errors 
    nErr_dqpsk_demod(m) = size(find([bin_data - estBit_demod]),2); % counting the number of errors

    % theoretical BER computation 
    a = sqrt(2*10.^(Eb_N0_dB(m)/10)*(1-sqrt(1/2)));
    b = sqrt(2*10.^(Eb_N0_dB(m)/10)*(1+sqrt(1/2)));
    
    k_bessel = 0:10; 
    temp = exp(-((a.^2+b.^2)/2)).*sum((a/b).^k_bessel.*besseli(k_bessel,a*b));
    theoryBer_dqpsk_demod(m) = temp - 0.5*besseli(0,a*b)*exp(-((a.^2+b.^2)/2));

end
% Simulation result for BER
simBer_dqpsk_demod = nErr_dqpsk_demod/N;

% generating pi/4-DQPSK modulated signal in time
mod_dqpsk = []; 
for g = 1: length(Phase_Gray)
    In_phase_signal = cos(2*pi*fc*t + pi*Phase_Gray(g)/4);
    Q_phase_signal = sin(2*pi*fc*t + pi*Phase_Gray(g)/4);
    mod_dqpsk = [mod_dqpsk (In_phase_signal+Q_phase_signal)];
end

% plotting pi/4-DQPSK modulated signal in time
dur = N*Tb;
tt = linspace(0,dur,length(mod_dqpsk));
figure,
subplot(211), stairs(bin_data,'LineWidth',1.5), grid on, ylim([0 1.25])
xlabel('\bf # of bits'), ylabel('\bf Amplitude'), title('Binary input')
subplot(212),plot(1000*tt, mod_dqpsk,'LineWidth',1.5),grid on, ylim([-1.5 1.5])
xlabel('\bf time (msec)'), ylabel('\bf Amplitude'), title('pi/4-DQPSK Modulated signal')

figure,
subplot(211), stairs(bin_data(1,1:20),'LineWidth',1.5), grid on, ylim([0 1.25])
xlabel('\bf # of bits'), ylabel('\bf Amplitude'), title('Binary input')
subplot(212),plot(1000*tt(1,1:1000), mod_dqpsk(1,1:1000),'LineWidth',1.5),grid on, ylim([-1.5 1.5])
xlabel('\bf time (msec)'), ylabel('\bf Amplitude'), title('pi/4-DQPSK Modulated signal')

%%
% plotting constellation of pi/4-DQPSK Modulation
figure,
plot(dqpsk_signal, 'bs-', 'Linewidth',2 , 'MarkerSize', 10, 'MarkerEdgeColor','r'), grid minor
axis([-1.25 1.25 -1.25 1.25])
xlabel('\bf In-phase'), ylabel('\bf Quadrature')
title('Constellation of pi/4-DQPSK Modulation')

% plotting constellation points of pi/4-DQPSK Modulation with/without AWGN
scatterplot(dqpsk_signal, [], [], 'r*'), grid on, title('\pi/4-DQPSK without AWGN')
scatterplot(y_dqpsk, [], [], 'r*'), grid on, title('\pi/4-DQPSK with AWGN')

%% Simulation of BER: computing over symbols 
% Set the modulation order to 4 to model DQPSK modulation.
M = 4;
% Apply DQPSK modulation to the input symbols. 
phase_DQPSK = pi/4;
y_sim_dqpsk = dpskmod(Bit_Gray, M, phase_DQPSK);
% Specify a constellation diagram
cd = comm.ConstellationDiagram('ShowTrajectory',true,'ShowReferenceConstellation',true);
cd(y_sim_dqpsk);
errorCalc = comm.ErrorRate;
totErr = 0;
numBits = 0;
mBER = zeros(length(SNR_dB),1);
for k = 1: length(SNR_dB)
    % Apply pi/4-DQPSK modulation to the input symbols
    y_sim_dqpsk = dpskmod(Bit_Gray, M, phase_DQPSK);
    rx_signal = awgn(y_sim_dqpsk,SNR_dB(1,k),'measured');
    decodmsg = dpskdemod(rx_signal,M,phase_DQPSK);      % Demodulate
    berVec = errorCalc(Bit_Gray,decodmsg);   % Calculate BER
    mBER(k, :) = berVec(1);
    totErr = totErr + berVec(2);
    numBits = numBits + berVec(3);
end

scatterplot(rx_signal), grid on

figure,
semilogy(Eb_N0_dB,mBER,'k*-', 'Linewidth',2), grid on
xlabel('\bf Eb/No [dB]')
ylabel('\bf Bit Error Rate')
legend('\pi/4-DQPSK simulation');

%% Simulation of BER: computing over bits
%
phase_DQPSK = pi/4;
txData = double(bin_data');
dpskmod = comm.DPSKModulator(M,phase_DQPSK,'BitInput',true);
dpskdemod = comm.DPSKDemodulator(M,phase_DQPSK,'BitOutput',true);
errorRate = comm.ErrorRate('ComputationDelay',1);
errorStats = zeros(length(Eb_N0_dB),3);
for k = 1:length(SNR_dB)
    channel = comm.AWGNChannel('EbNo',Eb_N0_dB(k),'BitsPerSymbol',k);
    modSig = dpskmod(txData);
    rxSig = channel(modSig);
    rxData = dpskdemod(rxSig);
    errorStats(k,:) = errorRate(txData,rxData);
end

scatterplot(rxSig), grid on

figure,
semilogy(Eb_N0_dB,errorStats(:,1),'k*-', 'Linewidth',2), grid on
xlabel('\bf Eb/No [dB]')
ylabel('\bf Bit Error Rate')
legend('\pi/4-DQPSK simulation');


%%

% plotting BER
figure,
semilogy(Eb_N0_dB,theoryBer_dqpsk_demod,'b-o', 'Linewidth',2);
hold on
semilogy(Eb_N0_dB,simBer_dqpsk_demod,'m*-','Linewidth',2);
% axis([0 14 10^-5 0.5])
grid on
legend('\pi/4-DQPSK Theoritical', '\pi/4-DQPSK simulation');
xlabel('\bf Eb/No [dB]')
ylabel('\bf Bit Error Rate')
title('probability of bit error curve for demodulation of \pi/4-DQPSK')

