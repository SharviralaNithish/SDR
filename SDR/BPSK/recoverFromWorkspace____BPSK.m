%save the data and assign it to a variable Hello.
alpha = Hello(1:10000,1,5);

prb = [0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1]';

prbdet = comm.PreambleDetector(prb,'Input','Bit');

idx = prbdet(alpha);

a = alpha((idx(4)+1):(idx(5)-16)).';

stem(a)

b = reshape(a, 8, (idx(4)+1):(idx(5)-16)/8).';

b = num2str(b);

%c = char(bin2dec(b));

%type(c);

fprintf(char(bin2dec(b)));