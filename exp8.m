clc; close all;clear all;
N = 1000; SNR_dB = 10; M = 4; loop_bw = 0.01; true_phase = pi/3;
tx_symbols = exp(1j * (2 * pi * (0:M-1) / M));
tx_data = randi([0 M-1], N, 1);
tx_signal = tx_symbols(tx_data + 1);
noise = (1/sqrt(2*10^(SNR_dB/10))) * (randn(N,1) + 1j*randn(N,1));
rx_signal = tx_signal .* exp(1j * true_phase) + noise;
est_phase_ml = angle(sum(conj(tx_signal) .* rx_signal));
est_phase_pll = zeros(N,1);
current_phase = 0;
for n = 1:N
    phase_err = angle(rx_signal(n) * exp(-1j * current_phase));
    current_phase = current_phase + loop_bw * phase_err;
    est_phase_pll(n) = current_phase;
end
corrected_rx_signal = rx_signal .* exp(-1j * est_phase_pll);
figure;
subplot(2,1,1); scatter(real(rx_signal), imag(rx_signal), 'filled');
title('Received Signal'); xlabel('In-Phase'); ylabel('Quadrature'); axis equal;
subplot(2,1,2); scatter(real(corrected_rx_signal), imag(corrected_rx_signal), 'filled');
title('Corrected Signal'); xlabel('In-Phase'); ylabel('Quadrature'); axis equal;
SNR_range = 0:2:20; phase_err_var = zeros(length(SNR_range),1);
for idx = 1:length(SNR_range)
    noise = (1/sqrt(2*10^(SNR_range(idx)/10))) * (randn(N,1) + 1j*randn(N,1));
    rx_signal = exp(1j * true_phase) + noise;
    phase_err_var(idx) = var(angle(rx_signal) - true_phase);
end
figure; plot(SNR_range, phase_err_var); title('Noise Effect on Phase Estimation');
xlabel('SNR (dB)'); ylabel('Phase Error Variance');
est_phase_dd = zeros(N,1); est_phase_ndd = zeros(N,1);
phase_dd = 0; phase_ndd = 0;
for n = 1:N
    noise = (1/sqrt(2*10^(SNR_dB/10))) * (randn + 1j*randn);
    rx_signal = tx_signal(n) * exp(1j * true_phase) + noise;
    detected_sym = exp(1j * round(angle(rx_signal) * M / (2*pi)) * 2*pi/M);
    phase_dd = phase_dd + loop_bw * angle(detected_sym * exp(-1j * phase_dd));
    phase_ndd = phase_ndd + loop_bw * angle(rx_signal * exp(-1j * phase_ndd));
    est_phase_dd(n) = phase_dd;
    est_phase_ndd(n) = phase_ndd;
end
figure;
subplot(2,1,1); plot(1:N, est_phase_dd); title('Decision-Directed Phase Estimate'); xlabel('Samples'); ylabel('Phase (radians)');
subplot(2,1,2); plot(1:N, est_phase_ndd); title('Non-Decision-Directed Phase Estimate'); xlabel('Samples'); ylabel('Phase (radians)');
