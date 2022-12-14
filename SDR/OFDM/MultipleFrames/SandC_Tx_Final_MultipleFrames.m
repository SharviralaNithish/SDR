 %% Maximum, 120 OFDM symbols(excluding Header) can be transmitted in one frame.
 clear all
 close all
 %NumFrames = 1;
 %% Build OFDM Modulator
 FFTLength = 64;
 NumGuardBandCarriers = [6; 5];
 NumDataCarriers = 48;
 CyclicPrefixLength = 16;
 PilotCarrierIndices = [12;26;40;54];
 NumOFDMSymInPreamble = 5;
 %NumBitsPerCharacter = 8;
 %% Convert message to bits

 Text = dec2bin('Hiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\nHiii Nithish. How are you? We are communicating from Mars.\n',8);
 numChars = size(Text,1);
 fprintf(char(bin2dec(Text)))

 % this line converts binary strings into vectors
 Text = Text - '0';
 % reshape into a vector of bits; 8 bits per character
 Text = reshape(Text.',8*numChars,1);
 %msgInBits = Text %repmat(randi([0 1], NumDataCarriers, 1),10, 1)
 %PayloadBits = msgInBits(:)
 % Calculate number of OFDM symbols without Header
 NumOFDMSymbolsWithOutHeader = ceil(length(Text)/NumDataCarriers);
 % Calculate number of bits padded in each frame
 NumPadBits = NumDataCarriers * NumOFDMSymbolsWithOutHeader - length(Text);
 % Get preamble for each frame
 Preamble = double(getOFDMPreambleAndPilot('Preamble',FFTLength, NumGuardBandCarriers));
 % Get pilot for each frame
 Pilots = double(getOFDMPreambleAndPilot('Pilot', 1));
 Text = [Text;randi([0,1],NumPadBits,1)];
 % make sure string is properly recovered
 % Rmessage = reshape(a(:), 8, numChars).'
 % Rmessage = num2str(Rmessage)
 % fprintf(char(bin2dec(Rmessage)))

 %% Build transmit packet
 %pre = [0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1]
 %message = [pre a]
 %% Build Header
 Header = dec2bin(numChars,32);
 Header = Header-'0';
 Header = Header(:);
 prb = [0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1]';
 Header = [prb;Header]
 NumFrames = ceil(length(Text)/4800);

% % message = zeros((length(Text)+length(Header)),1)
% %  if length(Text) > 4800 % Considering 100 OFDM symbols in each frame.
% %      message(1:48) = Header
% %      for k = 1:(NumFrames-1)
% %          if k == 1
% %              message(49:4848) = Text(1:4800)
% %          else
% %              message((4848+(4800*(k-2)+1)):(4848+(4800*(k-1)))) = Text((4800*(k-1)+1):(4800*(k)))
% %          end        
% %      end
% %      k = NumFrames
% %      message(((4848+(4800*(k-2)+1)):(4848+(4800*(k-2))+(NumOFDMSymbolsWithOutHeader-100*(NumFrames-1))*(48)))) = Text((4800*(k-1)+1):(4800*(k-1))+(NumOFDMSymbolsWithOutHeader-100*(NumFrames-1))*48)
% %  end
% % 
% %  if length(Text) < 4800
% %      message = [Header;Text]
% %  end
 message = [Header;Text];
 %numSymsTot = length(message)

 % BPSK modulator
 BPSKMod = comm.BPSKModulator;
 % OFDM modulator
 DataOFDMMod = comm.OFDMModulator(...
 'FFTLength' , FFTLength, ...
 'NumGuardBandCarriers', NumGuardBandCarriers, ...
 'InsertDCNull', true, ...
 'PilotInputPort', true, ...
 'PilotCarrierIndices', PilotCarrierIndices, ...
 'CyclicPrefixLength', CyclicPrefixLength, ...
 'NumSymbols', 1);

 %% Modulate Data
 symPostBPSK = BPSKMod.step(message);% Pad with preamble sequences instead so that they can be identified and removed. 
 % OFDM modulation for one frame
 symPostOFDM = zeros(80,NumOFDMSymbolsWithOutHeader+1);
 for i = 1:NumOFDMSymbolsWithOutHeader+1
     symPostEachOFDM = DataOFDMMod.step(symPostBPSK((48*(i-1)+1):(48*(i))), Pilots);
     symPostOFDM(:,i) = symPostEachOFDM;
 end 
 symPostOFDM = symPostOFDM(:);
 if NumFrames == 1
     y = [Preamble; symPostOFDM];
 else
     y = zeros(length(symPostOFDM)+length(Preamble)*NumFrames,1);
     for FrameNo = 1:NumFrames-1
         if FrameNo == 1
             y(1:8400) = [Preamble;symPostOFDM(1:8080)];
         else  
             y(8400+1+8320*(FrameNo-2):8400+8320*(FrameNo-1)) = [Preamble;symPostOFDM(8080+1+8000*(FrameNo-2):8080+8000*(FrameNo-1))];
         end
     end
     FrameNo = NumFrames;
     y(8400+1+8320*(FrameNo-2):(8400+(8320*(FrameNo-2))+320+(NumOFDMSymbolsWithOutHeader-100*(NumFrames-1))*(80))) = [Preamble;symPostOFDM(8080+1+8000*(FrameNo-2):8080+8000*(FrameNo-2)+(NumOFDMSymbolsWithOutHeader-100*(NumFrames-1))*80)];
 end
         
 % Repeat the frame
 %y = repmat([Preamble; symPostOFDM], NumFrames, 1)
 %y = repmat([randi([0,1],59,1);y],30,1) %For simulations.
 SamplesPerFrame = size(y,1);
 SampleTime = 1/SamplesPerFrame;    

% % Rdata = y(321:400)
% % DataOFDMMod = comm.OFDMModulator(...
% %  'FFTLength' , FFTLength, ...
% %  'NumGuardBandCarriers', NumGuardBandCarriers, ...
% %  'InsertDCNull', true, ...
% %  'PilotInputPort', true, ...
% %  'PilotCarrierIndices', PilotCarrierIndices, ...
% %  'CyclicPrefixLength', CyclicPrefixLength, ...
% %  'NumSymbols', 1)
% % odd = comm.OFDMDemodulator(DataOFDMMod)
% % [dataFreq,pilots] = odd(Rdata)
% % %Rmessage = dataFreq(:)
% % bpskdemodulator = comm.BPSKDemodulator 
% % symPostBPSKDemod = bpskdemodulator(dataFreq(:))
% % Rmessage = reshape(symPostBPSKDemod, 48, 1).'
% % Rmessage = num2str(Rmessage)
% % Rmessage = bin2dec(Rmessage)
% % %fprintf(char(bin2dec(Rmessage)))
% % 
% % Rdata = y(401:400+(Rmessage/6)*80)
% % DataOFDMMod = comm.OFDMModulator(...
% %  'FFTLength' , FFTLength, ...
% %  'NumGuardBandCarriers', NumGuardBandCarriers, ...
% %  'InsertDCNull', true, ...
% %  'PilotInputPort', true, ...
% %  'PilotCarrierIndices', PilotCarrierIndices, ...
% %  'CyclicPrefixLength', CyclicPrefixLength, ...
% %  'NumSymbols',Rmessage/6)
% % odd = comm.OFDMDemodulator(DataOFDMMod)
% % [dataFreq,pilots] = odd(Rdata)
% % %Rmessage = dataFreq(:)
% % bpskdemodulator = comm.BPSKDemodulator 
% % symPostBPSKDemod = bpskdemodulator(dataFreq(:))
% % Rmessage = reshape(symPostBPSKDemod, 8, Rmessage).'
% % Rmessage = num2str(Rmessage)
% % %Rmessage = bin2dec(Rmessage)
% % fprintf(char(bin2dec(Rmessage)))



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