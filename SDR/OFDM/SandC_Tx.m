 clear all
 close all
 NumFrames = 2;
 %% Build OFDM Modulator
 FFTLength = 64;
 NumGuardBandCarriers = [6; 5];
 NumDataCarriers = 48;
 CyclicPrefixLength = 16;
 PilotCarrierIndices = [12;26;40;54];
 NumOFDMSymInPreamble = 5;
 %NumBitsPerCharacter = 8;
 %% Convert message to bits

 Text = dec2bin('Hi ShNithish. How are you? We are communicating from Mars.\n',8)
 numChars = size(Text,1)
 fprintf(char(bin2dec(Text)))

 % this line converts binary strings into vectors
 Text = Text - '0'

 % reshape into a vector of bits; 8 bits per character
 Text = reshape(Text.',8*numChars,1)'

 % make sure string is properly recovered
 % Rmessage = reshape(a(:), 8, numChars).'
 % Rmessage = num2str(Rmessage)
 % fprintf(char(bin2dec(Rmessage)))

 %% Build transmit packet
 %pre = [0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1]
 %message = [pre a]
 message = Text
 numSymsTot = length(message)
 %%
 msgInBits = message %repmat(randi([0 1], NumDataCarriers, 1),10, 1)
 PayloadBits = msgInBits(:)
 % Calculate number of OFDM symbols per frame
 NumOFDMSymbols = ceil(length(PayloadBits)/NumDataCarriers)
 % Calculate number of bits padded in each frame
 NumPadBits = NumDataCarriers * NumOFDMSymbols - length(PayloadBits)
 % Get preamble for each frame
 Preamble = double(getOFDMPreambleAndPilot('Preamble',FFTLength, NumGuardBandCarriers))
 % Get pilot for each frame
 Pilots = double(getOFDMPreambleAndPilot('Pilot', NumOFDMSymbols))
 % BPSK modulator
 BPSKMod = comm.BPSKModulator
 % OFDM modulator
 DataOFDMMod = comm.OFDMModulator(...
 'FFTLength' , FFTLength, ...
 'NumGuardBandCarriers', NumGuardBandCarriers, ...
 'InsertDCNull', true, ...
 'PilotInputPort', true, ...
 'PilotCarrierIndices', PilotCarrierIndices, ...
 'CyclicPrefixLength', CyclicPrefixLength, ...
 'NumSymbols', NumOFDMSymbols)

 %% Modulate Data
 symPostBPSK = BPSKMod.step([PayloadBits; randi([0 1], NumPadBits, 1)])
 % OFDM modulation for one frame
 symPostOFDM = DataOFDMMod.step(reshape(symPostBPSK, ...
 NumDataCarriers, NumOFDMSymbols), Pilots)
 % Repeat the frame
 y = repmat([Preamble; symPostOFDM], NumFrames, 1)


Rdata = y(321:1120)
odd = comm.OFDMDemodulator(DataOFDMMod)
[dataFreq,pilots] = odd(Rdata)
%Rmessage = dataFreq(:)
bpskdemodulator = comm.BPSKDemodulator 
symPostBPSKDemod = bpskdemodulator(dataFreq(:))
Rmessage = reshape(symPostBPSKDemod, 8, numChars).'
Rmessage = num2str(Rmessage)
fprintf(char(bin2dec(Rmessage)))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig = getOFDMPreambleAndPilot(sigType, varargin)
% GETOFDMPREAMBLEANDPILOT: Return either the preamble OFDM symbols or the
% pilot signals for one frame transmission based on the 802.11a standard.

% Copyright 2014-2016 The MathWorks, Inc.

switch sigType
  case 'Preamble'
    [FFTLen, numGuardBandCarriers] = deal(varargin{:});
    
    % Create short preamble
    shortPreamble = ...
        [ 0    0  1+1i 0  0    0 -1-1i 0  0    0 ...                % [-27: -17]
          1+1i 0  0    0 -1-1i 0  0    0 -1-1i 0 0 0 1+1i 0 0 0 ... % [-16: -1 ]
          0    0  0    0 -1-1i 0  0    0 -1-1i 0 0 0 1+1i 0 0 0 ... % [ 0 :  15]
          1+1i 0  0    0  1+1i 0  0    0  1+1i 0 0 ].';             % [ 16:  27]

    % Create long Preamble
    longPreamble = complex(...
        [ 1,  1, -1, -1,  1,  1, -1,  1, -1,  1,  1,  1,...
          1,  1,  1, -1, -1,  1,  1, -1,  1, -1,  1,  1,  1,  1, 0,...
          1, -1, -1,  1,  1, -1,  1, -1,  1, -1, -1, -1, -1, -1,...
          1,  1, -1, -1,  1, -1,  1, -1,  1,  1,  1,  1].', 0);

    % OFDM preamble modulator
    preambleOFDMMod = comm.OFDMModulator(...
        'FFTLength' ,           FFTLen,...
        'NumGuardBandCarriers', numGuardBandCarriers,...
        'CyclicPrefixLength',   0);

    % Modulate and scale short preamble
    shortPreamblePostOFDM = sqrt(13/6)*preambleOFDMMod(shortPreamble);
    % Modulate long preamble
    longPreamblePostOFDM = preambleOFDMMod(longPreamble);

    % Preamble for one frame
    sig = [shortPreamblePostOFDM; ...
           shortPreamblePostOFDM; ...
           shortPreamblePostOFDM(1:end/2); ...
           longPreamblePostOFDM(end/2+1:end); ...
           longPreamblePostOFDM; ...
           longPreamblePostOFDM];
       
  case 'Pilot'
    numOFDMSym = varargin{1};
    
    % PN sequence for pilot generation
    obj.pPNSeq = comm.PNSequence(...
        'Polynomial',        [1 0 0 0 1 0 0 1],...
        'SamplesPerFrame',   numOFDMSym,...
        'InitialConditions', [1 1 1 1 1 1 1]);

    % Pilot signals for one frame
    pilotForOneCarrier  = obj.pPNSeq();
    sig = [repmat(1 - 2 * pilotForOneCarrier, 1, 3), ...
                      2 * pilotForOneCarrier - 1]';
end

end

% [EOF]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%