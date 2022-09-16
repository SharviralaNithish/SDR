%% generateMessage.m
%   Generates binary bit sequence for transmitting ASCII text with Ettus B200 USRP.
%
%   Cory J. Prust, Ph.D.
%   Last Modified:  7/18/2018

clear all;
close all;

%% Create message.
a = dec2bin('Hi Nithish. How are you? We are communicating from Mars.\n',8);
numChars = size(a,1);
fprintf(char(bin2dec(a)));

% this line converts binary strings into vectors
a = a - '0';

% reshape into a vector of bits; 8 bits per character
a = reshape(a.',8*numChars,1)';

% make sure string is properly recovered
b = reshape(a, 8, numChars).';
b = num2str(b);
fprintf(char(bin2dec(b)));

%% Build transmit packet
pre = [0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1];
message = [pre a];
numSymsTot = length(message);
