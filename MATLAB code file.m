clc; clear; close all;

%% Simulation Parameters
snrRange = 0:2:30;           % SNR range in dB
numBits = 1e6;               % Number of bits per SNR point
qamOrders = [16, 256];       % QAM modulation orders
codeRate = 1/2;              % Coding rate
maxIterations = 50;          % Decoding iterations

%% LDPC Configuration (DVB-S.2 Standard)
H = dvbs2ldpc(codeRate);     % Parity-check matrix
ldpcEncoder = ldpcEncoderConfig(H);
ldpcDecoder = ldpcDecoderConfig(H); % Critical fix

%% Turbo Code Configuration (3GPP Standard)
trellis = poly2trellis(4, [13 15], 13);
interleaver = randperm(ldpcEncoder.NumInformationBits);
turboEnc = comm.TurboEncoder('TrellisStructure', trellis, ...
                            'InterleaverIndices', interleaver);
turboDec = comm.TurboDecoder('TrellisStructure', trellis, ...
                            'InterleaverIndices', interleaver, ...
                            'NumIterations', maxIterations);

%% Initialize BER Storage
berUncoded = zeros(length(qamOrders), length(snrRange));
berLDPC = zeros(1, length(snrRange));
berTurbo = zeros(1, length(snrRange));

%% Main Simulation
for modIdx = 1:length(qamOrders)
    M = qamOrders(modIdx);
    for snrIdx = 1:length(snrRange)
        snr = snrRange(snrIdx);
        fprintf('Processing %d-QAM at %d dB...\n', M, snr);
        
        %% Uncoded QAM System
        txBits = randi([0 1], numBits, 1);
        txSym = qammod(txBits, M, 'InputType', 'bit', 'UnitAveragePower', true);
        rxSym = awgn(txSym, snr, 'measured');
        rxBits = qamdemod(rxSym, M, 'OutputType', 'bit', 'UnitAveragePower', true);
        [~, berUncoded(modIdx, snrIdx)] = biterr(txBits, rxBits);
        
        %% LDPC-Coded System (QAM-16 only)
        if M == 16
            % LDPC Encoding
            msg = randi([0 1], ldpcEncoder.NumInformationBits, 1);
            enc = ldpcEncode(msg, ldpcEncoder);
            
            % Ensure proper length for modulation
            if mod(length(enc), log2(M)) ~= 0
                enc = [enc; zeros(log2(M)-mod(length(enc),log2(M)), 1)];
            end
            
            % Modulation
            txSym = qammod(enc, M, 'InputType', 'bit', 'UnitAveragePower', true);
            
            % AWGN Channel
            rxSym = awgn(txSym, snr, 'measured');
            
            % LLR Demodulation
            noiseVar = 10^(-snr/10);
            llr = qamdemod(rxSym, M, 'OutputType', 'llr', ...
                          'UnitAveragePower', true, 'NoiseVariance', noiseVar);
            
            % LDPC Decoding
            dec = ldpcDecode(llr, ldpcDecoder, maxIterations);
            [~, berLDPC(snrIdx)] = biterr(msg, dec(1:ldpcEncoder.NumInformationBits));
            
            %% Turbo-Coded System
            encTurbo = turboEnc(msg);
            txSymTurbo = qammod(encTurbo, M, 'InputType', 'bit', 'UnitAveragePower', true);
            rxSymTurbo = awgn(txSymTurbo, snr, 'measured');
            llrTurbo = qamdemod(rxSymTurbo, M, 'OutputType', 'llr', ...
                               'UnitAveragePower', true, 'NoiseVariance', noiseVar);
            decTurbo = turboDec(llrTurbo);
            [~, berTurbo(snrIdx)] = biterr(msg, decTurbo);
        end
    end
end

%% Generate Plots
figure('Position', [100 100 1200 800], 'Name', 'Performance Analysis');

% QAM-16/256 BER Plot
subplot(2,2,1);
semilogy(snrRange, berUncoded(1,:), 'r-o', 'LineWidth', 2); hold on;
semilogy(snrRange, berUncoded(2,:), 'b-s', 'LineWidth', 2);
plot(25, 1e-6, 'kp', 'MarkerSize', 15, 'LineWidth', 2);
grid on; title('QAM-16/256 BER Performance');
xlabel('SNR (dB)'); ylabel('BER'); 
legend('16-QAM', '256-QAM', '1e^{-6} Reference', 'Location', 'best');
ylim([1e-7 1]);

% Error Reduction Plot
subplot(2,2,2);
reduction = berUncoded(1,:)./berLDPC;
bar(snrRange, reduction, 'FaceColor', [0.5 0.9 0.5]);
hold on; plot(snrRange, 1.4*ones(size(snrRange)), 'r--', 'LineWidth', 2);
title('LDPC Error Reduction (ITU-T G.975.1)');
xlabel('SNR (dB)'); ylabel('BER Reduction Factor');
grid on; ylim([1 2]);

% Code Comparison Plot
subplot(2,2,[3 4]);
semilogy(snrRange, berUncoded(1,:), 'k--', 'LineWidth', 2); hold on;
semilogy(snrRange, berLDPC, 'b-o', 'LineWidth', 2);
semilogy(snrRange, berTurbo, 'r-s', 'LineWidth', 2);
title('LDPC vs Turbo Codes Performance Comparison');
xlabel('SNR (dB)'); ylabel('BER'); grid on;
legend('Uncoded QAM-16', 'LDPC-Coded', 'Turbo-Coded', 'Location', 'best');
ylim([1e-7 1]);

%% Constellation Plot
figure('Name', 'QAM-16 Constellation at 25 dB SNR');
txBits = randi([0 1], 1e4, 1);
txSym = qammod(txBits, 16, 'InputType', 'bit', 'UnitAveragePower', true);
rxSym = awgn(txSym, 25, 'measured');
scatter(real(rxSym), imag(rxSym), 40, 'filled', 'MarkerFaceAlpha', 0.3);
title('QAM-16 Constellation at 25 dB SNR');
xlabel('In-Phase Component'); ylabel('Quadrature Component');
grid on; axis equal; xlim([-1.5 1.5]); ylim([-1.5 1.5]);
