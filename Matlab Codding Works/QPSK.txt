clc;
clear all;

N=200;
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



SNR = 14;
Tx_awgn = awgn(mod,SNR, 'measured');
Tx_const_awgn = awgn(const, SNR,'measured');
Rx = Tx_awgn;
Rx_const = Tx_const_awgn;

scatterplot(const, [], [], 'g*'); grid minor, title('Without AWGN');
scatterplot(Rx_const, [], [], 'r*'); grid minor, title('With AWGN');

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



N=200;
SNRdB= 5:1:25; % SNR for simulation
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
    
    sig_Rx = S + noise; % Recived signal
    
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


figure;
semilogy(SNRdB, BER_theo,'r-')  
hold on
semilogy(SNRdB, BER,'k*')                                 
xlabel('SNR[dB]')                                    
ylabel('Bit Error Rate');                                         
legend('Theoretical', 'Simulated');
title('Probability of Bit Error for QPSK Modulation');
grid on;
hold off;

