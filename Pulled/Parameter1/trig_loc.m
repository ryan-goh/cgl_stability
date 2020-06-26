%%%measure the trigger distance location as speed c is varied, testing for
%%%O(eps^(-1/2)) behavior

jf = 47;
ind = jf-6;

SOL = cell(ind,1);
SOLC = SOL;

for number=1:ind
    file_name = sprintf('speedSolCeq%d.dat',number*50);
    SOL{number} = load(file_name);
end

BIF = load('speedomegaPaper.dat');
C = BIF(:,1);
W = BIF(:,end);
CDAT = C(50:50:(ind)*50);
WDAT = W(50:50:(ind)*50);



    %%Concatenate two sides of solutions
for jj=1:ind    
    S = SOL{jj};
    M = length(S(:,1));
    SC = zeros(2*M,4);
    SC(:,1) = [S(:,1) - 1; S(:,1)];
    for ii = 2:4
        SC(:,ii) = [S(:,ii);S(:,ii+3)];
    end
    SOLC{jj} = SC;
    
end


LL = 500;

