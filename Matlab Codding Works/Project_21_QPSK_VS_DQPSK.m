clc;
clear all;

N=20000;
data1 = randi([0 1], 1, N); % Getting 1 and 0 randomly


%input information 
NZR = 2*data1-1; 
data = reshape(NZR,2,N/2);
fc = 10.^3;
Tb = 1/fc;
t = 0:(Tb/99):Tb;

mod = [];
Ac = 1;
In_phase = [];
Quadrature = [];

for (r=1:N/2)
    In = data(1,r)*Ac*cos(2*pi*fc*t);
    Qd = data(2,r)*Ac*sin(2*pi*fc*t);
    In_phase = [In_phase In];
    Quadrature = [Quadrature Qd];
    mod = [mod In+Qd];
end

M = 4;
const = [];
for (n =1:N/2)
    if data(1,n) == -1 && data(2,n) == -1 
        C = exp(j*pi/4);
    elseif data(1,n) == -1 && data(2,n) == 1
        C = exp(j*3*pi/4);
    elseif data(1,n) == 1 && data(2,n) == -1
        C = exp(j*5*pi/4);
    elseif data(1,n) == 1 && data(2,n) == 1
        C = exp(j*7*pi/4);     
    end
    const = [const C];
end


SNR = 20;
Tx_awgn = awgn(mod,SNR, 'measured');
Tx_const_awgn = awgn(const, SNR,'measured');
Rx = Tx_awgn;
Rx_const = Tx_const_awgn;

scatterplot(const, [], [], 'g*'); grid minor, title('QPSK Without AWGN');
scatterplot(Rx_const, [], [], 'r*'); grid minor, title('QPSK With AWGN');

figure('Name','QPSK Modulation','NumberTitle','off');
subplot(2,1,1);
stairs(NZR);    
grid minor; 
xlim([0,N]); 
ylim([0,1.1]);
title('Message Signal');
subplot(2,1,2);
plot(mod); 
grid minor; 
title('QPSK Modulated Signal'); 
ylim([-2,2]);


SNRdB= 0:1:15; % SNR for simulation
SNRlin=10.^(SNRdB/10);
BER = zeros(1,length(SNRlin));% simulated BER
b1 = rand(1,N) > 0.5;
b2 = rand(1,N) > 0.5;
% QPSK symbol mapping
I = (2*b1) - 1;
Q = (2*b2) - 1;
S = I + 1j*Q; 
N0 = 1./SNRlin; % Variance
for k = 1:length(SNRdB)
    
    noise = sqrt(N0(k)/2)*(randn(1,N) + 1j*randn(1,N)); % AWGN noise
    
    sig_Rx = S + noise; % Received signal
    
    % For BER calculation
    sig_I = real(sig_Rx); % I component
    sig_Q = imag(sig_Rx); % Q component
    
    bld_I = sig_I > 0; % I decision 
    bld_Q = sig_Q > 0; % Q decision
    
    b1_error = (bld_I ~= b1); % Inphase bit error
    b2_error = (bld_Q ~= b2); % Quadrature bit error
    
    Error_bit = sum(b1_error) + sum(b2_error); % Total bit error
    BER(k) = sum(Error_bit)/(2*N); % Simulated BER
    
end
BER_theo = 2*qfunc(sqrt(2*SNRlin)); % Theoretical BER 

%%%%%%%%%%%%%%%%% DQPSK MODULATION &&&&&&&&&&&&&&&&&


M = 4; % QPSK constellation
k = log2(M); % number of bits per symbol

Eb_N0_dB = 0:1:15; % multiple Eb/N0 values
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

% plotting constellation of pi/4-DQPSK Modulation

figure,
plot(dqpsk_signal, 'bs-', 'Linewidth',2 , 'MarkerSize', 10, 'MarkerEdgeColor','r'), grid minor
axis([-1.25 1.25 -1.25 1.25])
xlabel('\bf In-phase'), ylabel('\bf Quadrature')
title('Constellation of pi/4-DQPSK Modulation')

% plotting constellation points of pi/4-DQPSK Modulation with/without AWGN

scatterplot(dqpsk_signal, [], [], 'r*'), grid on, title('\pi/4-DQPSK without AWGN')
scatterplot(y_dqpsk, [], [], 'r*'), grid on, title('\pi/4-DQPSK with AWGN')




figure;
semilogy(Eb_N0_dB,simBer_dqpsk_demod,'m*','Linewidth',2); %burasý
grid on
hold on
semilogy(SNRdB, BER,'g*','Linewidth',2);
xlabel('SNR[dB]');                                    
ylabel('Bit Error Rate');  
title('Probability of Bit Error for QPSK and DQPSK');
legend('DQPSK Simulation', 'QPSK Simulation');