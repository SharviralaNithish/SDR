for i = 1:size(Hello,3)

alpha = Hello(1:15000,1,i);

prb = [0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1]';

prbdet = comm.PreambleDetector(prb,'Input','Bit');

idx = prbdet(alpha);

if length(idx)>10


    a = alpha((idx(4)+1):(idx(5)-16)).';

    stem(a)

    b = reshape(a, 8,numChars).';

    b = num2str(b);

    fprintf(char(bin2dec(b)));
break

end
end
