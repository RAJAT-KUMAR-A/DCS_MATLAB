clc; close all;clear all;
num_bits = 20; samples_per_bit = 120; num_carriers = 6; samples = [10, 20, 30, 40, 60, 120];
disp('Enter your bit sequence:');
bit_sequence = str2num(input('', 's'));
if length(bit_sequence) ~= num_bits, error('Bit length mismatch!'); end
input_signal = repelem(2*bit_sequence - 1, samples_per_bit);
carrier_signal = cos(linspace(0, 2*pi*num_bits, samples_per_bit*num_bits));
bpsk_mod_signal = input_signal .* carrier_signal;
carriers = cell(1, num_carriers);
for i = 1:num_carriers
    t = linspace(0, 2*pi, samples(i) + 1); t(end) = []; 
    carriers{i} = repmat(cos(t), 1, ceil(samples_per_bit / length(t)));
    carriers{i} = carriers{i}(1:samples_per_bit);
end
spread_signal = [];
for i = 1:num_bits
    carrier_idx = randi([1, num_carriers]);
    spread_signal = [spread_signal carriers{carrier_idx}];
end
freq_hopped_sig = bpsk_mod_signal .* spread_signal;
bpsk_demodulated = freq_hopped_sig ./ spread_signal;
original_BPSK_signal = bpsk_demodulated ./ carrier_signal;
figure(1);
subplot(4,1,1); plot(input_signal); title('Original Bit Sequence');
subplot(4,1,2); plot(bpsk_mod_signal); title('BPSK Modulated Signal');
subplot(4,1,3); plot(spread_signal); title('Spread Signal with 6 Frequencies');
subplot(4,1,4); plot(freq_hopped_sig); title('Frequency Hopped Spread Signal');
figure(2);
subplot(2,1,1); plot(bpsk_demodulated); title('Demodulated BPSK Signal');
subplot(2,1,2); plot(original_BPSK_signal); title('Transmitted Original Bit Sequence');


