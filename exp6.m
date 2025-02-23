L        = 4;         
rollOff  = 0.5;      
rcDelay  = 10;  
htx = rcosdesign(rollOff, 6, 4);
hrx  = conj(fliplr(htx));
p = conv(htx,hrx);
M = 2; 
data = zeros(1, 2*rcDelay);
data(1:2:end) = 1;
txSym = real(pammod(data, M));
txUpSequence = upsample(txSym, L);
txSequence = filter(htx, 1, txUpSequence);
timeOffset = 1; 
rxDelayed = [zeros(1, timeOffset), txSequence(1:end-timeOffset)];
mfOutput = filter(hrx, 1, rxDelayed);
rxSym = downsample(mfOutput, L);
selectedSamples = upsample(rxSym, L);
selectedSamples(selectedSamples == 0) = NaN;
figure
plot(complex(rxSym(rcDelay+1:end)), 'o')
grid on
xlim([-1.5 1.5])
title('Rx Scatterplot')
xlabel('In-phase (I)')
ylabel('Quadrature (Q)')
figure
stem(rxSym)
title('Symbol Sequence with delay')
xlabel('Symbol Index')
ylabel('Amplitude')
rxSym = downsample(mfOutput, L, timeOffset);
selectedSamples = upsample(rxSym, L);
selectedSamples(selectedSamples == 0) = NaN;
figure
plot(complex(rxSym(rcDelay+1:end)), 'o')
grid on
xlim([-1.5 1.5])
title('Rx Scatterplot')
xlabel('In-phase (I)')
ylabel('Quadrature (Q)')
figure
stem(rxSym)
title('Symbol Sequence without delay')
xlabel('Symbol Index')
ylabel('Amplitude')

