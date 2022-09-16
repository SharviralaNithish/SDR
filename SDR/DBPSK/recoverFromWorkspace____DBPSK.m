alpha = zeros(10000,1);

for i = 1:9999
    if abs(wrapTo2Pi(angle(Hello(i,1,5)))-wrapTo2Pi(angle(Hello(i+1,1,5))))<1.57
       alpha(i) = 0;
    else
       alpha(i) = 1;
    end    
end     

prb = [0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1]';

prbdet = comm.PreambleDetector(prb,'Input','Bit');

idx = prbdet(alpha);

a = alpha((idx(4)+1):(idx(5)-16)).';

stem(a)

b = reshape(a, 8, numChars).';

b = num2str(b);

fprintf(char(bin2dec(b)));