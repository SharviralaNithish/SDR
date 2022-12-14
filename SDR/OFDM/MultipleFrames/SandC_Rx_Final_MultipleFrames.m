%clear all
%close all
FFTLength = 64;
NumGuardBandCarriers = [6; 5];
PilotCarrierIndices = [12;26;40;54];
CyclicPrefixLength = 16;
prb = [0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1]';
%NumOFDMSymbols = 11
%numChars = 66
Preamble = double(getOFDMPreambleAndPilot('Preamble',FFTLength, NumGuardBandCarriers));

%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%% %% %%
% Schmidl and Cox: Coarse Packet Detection 

 L = 16; % Short sync field length
 m = L; % Distance between fields
 N = 100000; % Autocorrelation samples
 M = zeros(N,1);

 for k1=1:N
 P = (conj(alpha(k1:k1+m-1)).')*alpha(k1+L:k1+m+L-1);
 dr = abs(alpha(k1+L:k1+m+L-1));
 R = dr'*dr;
 M(k1) = abs(P)/(R);
 end

 stem(M)
 grid on;
 xlabel('k')
 ylabel('M')
 legend('Autocorrelation')
 title('Packet detection')
%% Estimating the positions of packets.
 count = 0;
 target = 0;
 Est_Indices = zeros(30,1);
 for round=1:99998
     if ((0.99<M(round))&&(M(round)<1.01))&&((0.99<M(round+1))&&(M(round+1)<1.01))&&((0.99<M(round+2))&&(M(round+2)<1.01))
         count = count+1;
         if count>100
             target = target + 1;
             Est_Indices(target) = round;
             count = 0;
         end
     end
 end
%% Extracting the Packet
for i = 1:30
    r = alpha(Est_Indices(i)-100:Est_Indices(i+1));
    % Determine STO
    %%Estimate fine sample offset
    LSTF = Preamble(1:160);
    LLTF = Preamble(161:161+160-1);
    symLen = 80; % FFTLength + CPLen---(refering to what?)
    known = LLTF(1:symLen,1);
    % Filter
    coeff = conj(flipud(known));
    c_filt = abs(filter(coeff, 1, r)).^2;
    % Correlation
    crosscorrelation = abs(xcorr(r,known)).^2;
    padding = length(r)-symLen+1; % Remove padding added to known sequence
    c_cor = crosscorrelation(padding:end);
    % Plot comparison
    stem(c_filt)
    hold on; stem(c_cor,'r*')
    [v1,peak1] = max(c_cor);
    c_cor(peak1) = 0;
    [v2,peak2] = max(c_cor);
    v = min([v1 v2]); % Does not make any sense as v2 will always be less than v1.
    plot([peak1 peak2],[v v]*0.9,'k','LineWidth',2);hold off;
    text(peak1+100,v*0.9,'64 Sample Gap')
    grid on;
    xlabel('Samples');ylabel('Magnitude')
    % Get numerical offset
    if abs(peak2-peak1)==FFTLength
        % Adjust position from start of LLTF
        minpeak = min([peak1 peak2]);
        LLTF_Start_est = minpeak-symLen;
        LSTF_Start_est = LLTF_Start_est - length(LSTF);
        % Display
        disp(LSTF_Start_est)
    end
     
    %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
    % Equalization Example
    r = r(LSTF_Start_est+1:end); % Removing offset
    %% Channel estimation
    preambleOFDMMod = comm.OFDMModulator(...
        'FFTLength' ,           FFTLength,...
        'NumGuardBandCarriers', NumGuardBandCarriers,...
        'CyclicPrefixLength',   0,'NumSymbols', 2,'InsertDCNull', true);
    od = comm.OFDMDemodulator(preambleOFDMMod);
    od.PilotOutputPort = true;     %%(???)
    % OFDM Demodulate LLTF
    LLTF = Preamble(161:161+160-1);
    rLLTF = r(161+32:161+160-1);
    [rLLTFFreq,rp] = od(rLLTF);
    [LLTFFreq,p] = od(LLTF(33:end));% remove CP
    % Estimate channel
    ls = rLLTFFreq./LLTFFreq; % Least-square estimate
    chanEst = mean(ls,2); % Average over both symbols
    CSI = real(chanEst.*conj(chanEst));% (CSI-Channel state information)
    ls = rp./p; % Least-square estimate
    chanEstPilots = mean(ls,2); % Average over both symbols
    CSIPilots = real(chanEstPilots.*conj(chanEstPilots));  % (CSI-Channel state information)
    %% Demodulating Header
    data = r(321:400);
    DataOFDMMod = comm.OFDMModulator(...
     'FFTLength' , FFTLength, ...
     'NumGuardBandCarriers', NumGuardBandCarriers, ...
     'InsertDCNull', true, ...
     'PilotInputPort', true, ...
     'PilotCarrierIndices', PilotCarrierIndices, ...
     'CyclicPrefixLength', CyclicPrefixLength, ...
     'NumSymbols', 1);
    odd = comm.OFDMDemodulator(DataOFDMMod);
    [dataFreq,pilots] = odd(data);
    % Apply LLTF's estimate to data symbols and data pilots
    postLLTFEqData = bsxfun(@times, dataFreq, conj(chanEst(:))./CSI(:));
    postLLTFEqPilots = ...
        bsxfun(@times, pilots, conj(chanEstPilots(:))./CSIPilots(:));
    % Visualization objects
    tt1 = comm.ConstellationDiagram;
    tt2 = comm.ConstellationDiagram;
    tt2.Position = tt2.Position + [500 0 0 0];
    % Estimate remaining offsets with pilots
    correctedSymbols = zeros(size(postLLTFEqData));
    for symbol = 1:size(postLLTFEqData,2)
        % Estimate rotation across pilots
        p = postLLTFEqPilots(:,symbol);
        e = conj(mean(p.*conj(Pilots(:,symbol))));
        % Equalize
        sig = postLLTFEqData(:,symbol).*e;
        correctedSymbols(:,symbol) = sig;
        % Visualize
        tt1(sig);tt2(postLLTFEqData(:,symbol));pause(0.1)
    end
    bpskdemodulator = comm.BPSKDemodulator; 
    symPostBPSKDemod = bpskdemodulator(double(correctedSymbols(:)));

    if symPostBPSKDemod(1:16) == prb
        break
    end
end   
Header = (symPostBPSKDemod(17:48)).';
Header = num2str(Header);
Header = bin2dec(Header);
NumOFDMSymbolsAfterHeader = ceil(Header/6);
NumFrames = ceil(NumOFDMSymbolsAfterHeader/100);
%% Extracting the entire packet                                
data = alpha(Est_Indices(i)-100+LSTF_Start_est:Est_Indices(i)-100+LSTF_Start_est-1+(length(Preamble)*NumFrames)+(NumOFDMSymbolsAfterHeader+1)*80); 

if NumFrames == 1
    correctedSymbols = zeros(48,NumOFDMSymbolsAfterHeader);
    for j = 1:NumOFDMSymbolsAfterHeader
        [dataFreq,pilots] = odd(data((400+80*(j-1)+1):400+(80*(j))));
        % Apply LLTF's estimate to data symbols and data pilots
        postLLTFEqData = bsxfun(@times, dataFreq, conj(chanEst(:))./CSI(:));
        postLLTFEqPilots = ...
              bsxfun(@times, pilots, conj(chanEstPilots(:))./CSIPilots(:));
        p = postLLTFEqPilots(:);
        e = conj(mean(p.*conj(Pilots(:))));
        % Equalize
        sig = postLLTFEqData(:).*e;
        correctedSymbols(:,j) = sig;
        % Visualize
        %%%tt1(sig);tt2(postLLTFEqData(:,symbol));pause(0.1)
    end
else
    correctedSymbols = zeros(48,NumOFDMSymbolsAfterHeader);
    % Demodulating first frame
    for j = 1:100
        [dataFreq,pilots] = odd(data((400+80*(j-1)+1):400+(80*(j))));
        % Apply LLTF's estimate to data symbols and data pilots
        postLLTFEqData = bsxfun(@times, dataFreq, conj(chanEst(:))./CSI(:));
        postLLTFEqPilots = ...
              bsxfun(@times, pilots, conj(chanEstPilots(:))./CSIPilots(:));
        p = postLLTFEqPilots(:);
        e = conj(mean(p.*conj(Pilots(:))));
        % Equalize
        sig = postLLTFEqData(:).*e;
        correctedSymbols(:,j) = sig;
    end
    % Estimating subsequent channels and demodulating.
    for FrameNo = 2:NumFrames
        if FrameNo == NumFrames
            break
        end  
        %% Channel estimation
        preambleOFDMMod = comm.OFDMModulator(...
            'FFTLength' ,           FFTLength,...
            'NumGuardBandCarriers', NumGuardBandCarriers,...
            'CyclicPrefixLength',   0,'NumSymbols', 2,'InsertDCNull', true);
        od = comm.OFDMDemodulator(preambleOFDMMod);
        od.PilotOutputPort = true;     %%(???)
        % OFDM Demodulate LLTF
        LLTF = Preamble(161:161+160-1);
        rLLTF = data(8400+8320*(FrameNo-2)+161+32:8400+8320*(FrameNo-2)+161+160-1);
        [rLLTFFreq,rp] = od(rLLTF);
        [LLTFFreq,p] = od(LLTF(33:end));% remove CP
        % Estimate channel
        ls = rLLTFFreq./LLTFFreq; % Least-square estimate
        chanEst = mean(ls,2); % Average over both symbols
        CSI = real(chanEst.*conj(chanEst));% (CSI-Channel state information)
        ls = rp./p; % Least-square estimate
        chanEstPilots = mean(ls,2); % Average over both symbols
        CSIPilots = real(chanEstPilots.*conj(chanEstPilots));  % (CSI-Channel state information)
        for j = (FrameNo-1)*100+1:(FrameNo-1)*100+100
            [dataFreq,pilots] = odd(data((8400+320*(FrameNo-1)+80*(j-1-100)+1):8400+320*(FrameNo-1)+(80*(j-100))));
            % Apply LLTF's estimate to data symbols and data pilots
            postLLTFEqData = bsxfun(@times, dataFreq, conj(chanEst(:))./CSI(:));
            postLLTFEqPilots = ...
                  bsxfun(@times, pilots, conj(chanEstPilots(:))./CSIPilots(:));
            p = postLLTFEqPilots(:);
            e = conj(mean(p.*conj(Pilots(:))));
            % Equalize
            sig = postLLTFEqData(:).*e;
            correctedSymbols(:,j) = sig;
        end
    end     
        %% Channel estimation for last frame
        FrameNo = NumFrames;
        preambleOFDMMod = comm.OFDMModulator(...
            'FFTLength' ,           FFTLength,...
            'NumGuardBandCarriers', NumGuardBandCarriers,...
            'CyclicPrefixLength',   0,'NumSymbols', 2,'InsertDCNull', true);
        od = comm.OFDMDemodulator(preambleOFDMMod);
        od.PilotOutputPort = true;     %%(???)
        % OFDM Demodulate LLTF
        LLTF = Preamble(161:161+160-1);
        rLLTF = data(8400+8320*(FrameNo-2)+161+32:8400+8320*(FrameNo-2)+161+160-1);
        [rLLTFFreq,rp] = od(rLLTF);
        [LLTFFreq,p] = od(LLTF(33:end));% remove CP
        % Estimate channel
        ls = rLLTFFreq./LLTFFreq; % Least-square estimate
        chanEst = mean(ls,2); % Average over both symbols
        CSI = real(chanEst.*conj(chanEst));% (CSI-Channel state information)
        ls = rp./p; % Least-square estimate
        chanEstPilots = mean(ls,2); % Average over both symbols
        CSIPilots = real(chanEstPilots.*conj(chanEstPilots));  % (CSI-Channel state information)
        for j = (FrameNo-1)*100+1:NumOFDMSymbolsAfterHeader
            [dataFreq,pilots] = odd(data((8400+320*(FrameNo-1)+80*(j-1-100)+1):8400+320*(FrameNo-1)+(80*(j-100))));
            % Apply LLTF's estimate to data symbols and data pilots
            postLLTFEqData = bsxfun(@times, dataFreq, conj(chanEst(:))./CSI(:));
            postLLTFEqPilots = ...
                  bsxfun(@times, pilots, conj(chanEstPilots(:))./CSIPilots(:));
            p = postLLTFEqPilots(:);
            e = conj(mean(p.*conj(Pilots(:))));
            % Equalize
            sig = postLLTFEqData(:).*e;
            correctedSymbols(:,j) = sig;
        end    
end 

bpskdemodulator = comm.BPSKDemodulator; 
symPostBPSKDemod = bpskdemodulator(double(correctedSymbols(:)));
UsefulData = symPostBPSKDemod(1:Header*8);
Header = reshape(UsefulData, 8, Header).';
Header = num2str(Header);
fprintf(char(bin2dec(Header)))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


