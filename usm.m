function [u,kf,kb,y]=usm(Seq,AlphaB,L,S,seed)

%USM Universal sequence maps
%Syntax: [u,kf,kb,y]=usm(Seq,AlphaB,L,S,seed)
%Description: Seq is the sequence and AlphB is the alphabet. If the
%             alphabet is not provided, or is empty, then the unique sequences
%             in Seq will be used. Seq can also be the fileneme where the
%             sequence is kept. In that case fasta format will be assumed.
%
%             u is a matrix with the USM coordinates, e.g. if u has
%             size n x m then the first n/2 rows contain the forward USM
%             coordinates (which is to say Chaos Game representatin, CGR, with a
%             seed that is different from 1/2) and the seccond n/2 rows
%             contain the backward coordinates (which is the same as CGR on
%             the reverse direction). Note that the steady state USM
%             solution to seeding causes the first and last USM coordinates
%             of the backward and forward iteration to be the same. Note
%             that the first forward coordinate is the first to be iterated
%             but first backward coordinate is the last to be iterated,
%             that is, the coordinates converge to the same values, this
%             would not happen with fixed value seeds (see Chaos Game 
%             Representation, CGR, of the same sequence at the very end 
%             of this help under "Misc Notes").
%
%             >> usm('ABBBAAB')
%             ans =
%             0.4471    0.7236    0.8618    0.9309    0.4654    0.2327    0.6164
%             0.4471    0.8943    0.7885    0.5770    0.1541    0.3082    0.6164
%
%             If a seccond output argument is used, kf, the fractal kernel
%             representation of density will be produced for the forward
%             USM coordinates. In that case this function expects an additional
%             two input arguments, L and S, for the memory length and
%             smoothing settings of the kernel. The default values are S=1
%             and L=2;
%
%             The third output argumnt, Kb, is just like Kf but for the
%             backward USM coordinates.
%
%             Note that the kernel densities are calculated to have an
%             integral that is equal to the number of units in the
%             sequence. For example, for the 8 unit sequence example:
%
%             >> [u,kf]=usm('ABBBAACB');sum(y.kf(:))/prod(size(y.kf))
%             ans =
%             8
%
%             A fourth output argument is allowed, y, which will be a
%             structured variable with all the information included as input
%             and output arguments. This thrid variable is in the format
%             that usm_plot expects for graphic displaying.
%
%             For example, the calculation of the frequency of all tripples
%             in a sequence from an alphabet with 5 units, which causes USM
%             space to have ceil(log2(5))=3, can be calculated by using
%             a large value for the smoothing parameter, S, say 10^10:
%
%             >> [u,kf,kb,y]=usm('ABBBAACDEB','',3,10^10);y
%
%             y =
%                Seq: {'ABBBAACDEB'}
%             AlphaB: 'ABCDE'
%                USM: {[6x10 double]}
%                 kf: [8x8x8 double]
%                 kb: [8x8x8 double]
%
%             the structured variable y retains all the information.
%
%             There is one last, the 5th, input argument, the seed for
%             the USM iteration. If not provided, the looping
%             implementation that turns USM coordinates into steady state
%             solutions will be used. Allowing this argument serves the
%             purpose of enabling the use of this function to implement the
%             standard CGR - by providing 1/2 as the seed - or any other
%             fixed value seeding variation such as the use of a random
%             starting point also used in published reports.
%
%             Misc Notes:
%
%             submitting empty arguments will invoque their default values.
%             For example, CGR of a sequence can be 
%             >> u=usm('ABBBAAB','',[],[],1/2)
%             u =
%             0.2500    0.6250    0.8125    0.9063    0.4531    0.2266    0.6133
%             0.4492    0.8984    0.7969    0.5938    0.1875    0.3750    0.7500
%             >> cgr_forward=u(1,:)
%             cgr_forward =
%             0.2500    0.6250    0.8125    0.9063    0.4531    0.2266    0.6133
%             
%             If a structured variable is submitted as firts intput
%             argument, Seq, it is assumed that it is a variable of the
%             same type as the forth output argument y, carrying th eUSM
%             coordinates already. This can come handy to generate
%             different densities, with different values of L and S for the
%             same sequence whitout having to recalculate de USM
%             coordinates, which stay the same. For example:
%             >> [u,kf,kb,y]=usm('ABBBAACDEB','',0,1);y   <--- stores USM values in y
%             y = 
%                Seq: {'ABBBAACDEB'}
%             AlphaB: 'ABCDE'
%                USM: {[6x10 double]}
%                 kf: 10
%                 kb: 10
%             kParms: [1x1 struct]
%             >> [u,kf,kb,y3]=usm(y,'',3,1);y3     <--- the USM coordinates are retrieved from y
%             y3 = 
%                Seq: {'ABBBAACDEB'}
%             AlphaB: 'ABCDE'
%                USM: {[6x10 double]}
%                 kf: [8x8x8 double]               <--- forward density for L=3
%                 kb: [8x8x8 double]
%             kParms: [1x1 struct]
%             >> [u,kf,kb,y4]=usm(y,'',4,1);y4     <--- the USM coordinates are retrieved from y
%             y4 = 
%                Seq: {'ABBBAACDEB'}
%             AlphaB: 'ABCDE'
%                USM: {[6x10 double]}
%                 kf: [16x16x16 double]             <--- forward density for L=4
%                 kb: [16x16x16 double]
%             kParms: [1x1 struct]
%
%             this is of course a silly example because the sequence is so
%             small that it is not worth the trouble of not calculating the
%             USM coordinates each time. It could however come handy for
%             very large sequences.
%
%Jonas Almeida 28 February 2006

%1. Check if Seq is a sequence or a text file with the sequence
if ischar(Seq)
    if exist(Seq)==2 %then sequence is to be read from fasta formated text file
        i=0;fid=fopen(Seq,'r');
        while ~feof(fid)
            line=fgetl(fid);
            if line(1)=='>'
                i=i+1;y.Header{i}=line(2:end);y.Seq{i}='';
            else
                y.Seq{i}=[y.Seq{i},line];
            end
        end
        fclose(fid);
    else
        y.Seq={Seq};
    end


    %2. Calculate usm coordinates

    %2.1. get the alphabet
    if nargin<2;AlphaB='';end
    if isempty(AlphaB); %then extract it from symbol usage by the sequences
        y.AlphaB=sort(unique([y.Seq{:}]));
    else
        y.AlphaB=AlphaB;
    end
    %2.2. generate binary sequence
    if nargin<5;seed='loop';end
    for i=1:length(y.Seq)
        B=Seq2Bin(y.Seq{i},y.AlphaB);
        %y.Bin{i}=B; %uncomment if you want to clutter y with the binary coordinates
        u_f=USM_CGR(B,seed); %forward iteration
        u_b=USM_CGR(B(:,end:-1:1),seed);u_b=u_b(:,end:-1:1); %backward iteration
        y.USM{i}=[u_f;u_b];
    end

    
elseif isstruct(Seq) % the structured variable is being submitted aleady
    y=Seq;
else
    error('argument format not recognized - it can only be char or struct')
end
if length(y.USM)==1
    u=y.USM{1};
else
    for i=1:length(y.USM)
        u{i}=y.USM{i};
    end
end
%3. Calculate Kernel densities
%Get Kernel parms

if nargin<3;L=[];end;if isempty(L);L=2;end
if nargin<4;S=[];end;if isempty(S);S=1;end

if nargout>1
    
    y.kf=usm_kernel(u(1:end/2,:),L,S);%forward densities
    kf=y.kf;
    if nargout>2
        y.kb=usm_kernel(u(end/2+1:end,:),L,S);%backward densities
        kb=y.kb;
        y.kParms.L=L;
        y.kParms.S=S;
    end
end

% ----------- NESTED FUNCTION ----------------
function B=Seq2Bin(S,A)

%Seq2Bin converts sequence in compact USM binary coordinates
%Syntax: B=Seq2Bin(S,A)
%Description: S is the sequence, A is teh alphabet and B is the binnary
%             matrix of USM coordinates.
%
%Jonas Almeida, 2 March 2006

%1. produce sparse USM coordinates
n=length(A);
for i=1:n
    B(i,:)=(S==A(i));
end
%2. convert sparce into compact coordinates
B=USM_compact(B);

% ----------- NESTED FUNCTION ----------------
function M=USM_compact(M)

%USM_COMPACT compacts hi matrix from sparce to compact USM
%Syntax: M=USM_compact(M)
%Description:
%  compacts a binary hit matrix, M. For example, compacting a
%  unidirectiopnal USM of a nucleotide sequence will produce the original
%  CGR hit matrix.
%
%Jonas Almeida, almeidaj@musc.edu, Nov 2004

[n,m]=size(M);
K=repmat([1:n]',1,m);
M=M.*K;
M=sum(M);
B=dec2bin(0:n-1)';
M=B(:,M)=='1';

% ----------- NESTED FUNCTION ----------------
function M=USM_CGR(M,seed)

%USM_CGR Chaos Game Representation (CGR) iteration for a Hit matrix
% Syntax: y=USM_CGR(M,seed)
% Description: using a binary hit matrix, M, this function iterates the
% forward CGR coordinates. By doing M=M(:,end:-1:1) and reversing the y
% too: y=y(:,end:-1:1),
% the optional seed allows:
%    seed='random'
%    seed='loop'     (default choice, uses the last coordinates as seed)
%    seed=0.5        (or anyother value, it assigns that value to all
%                    coordinates.A seed vector instead of a seed scalar can
%                    also be used)
%
%Jonas Almeida, almeidaj@musc.edu, Nov 2004


if nargin<2;seed='loop';end

M=double(M);
[n,m]=size(M);%y=zeros(n,m);

% All of this is worrying about the seed
if isnumeric(seed)
    CGRo=seed;
    if length(CGRo)==1;CGRo=ones(n,1).*CGRo;end

elseif ischar(seed)
    switch seed
        case 'random'
            CGRo=rand(n,1);
        case 'loop'
            %CGRo=USM_CGR(M(:,end-min([32,m-1]):end),0.5); % do CGR of the last 32 symbols or as many as available if less
            %CGRo=CGRo(:,end);
            %If the sequence is shorter than 2000 units then extend it with
            %loops to find steady state solution with the numerical
            %resolution of the processor.
            if m<10000
                if m==1; % If there is only one unit
                    Mo=repmat(M,1,2); %then build seed from dimer --> try it, you'll see why
                else
                    Mo=M;
                end
                MM=[];
                for i=1:ceil(10000/m)
                    MM=[MM,Mo,Mo(:,end-1:-1:2)];
                end
            else
                MM=M;
            end
            CGRo=USM_CGR(MM(:,9999:-1:2),0.5); % do CGR of the last 32 symbols or as many as available if less
            CGRo=CGRo(:,end);

        otherwise
            error('CGR seeding option not defined')
    end
else
    error(['seed variable needs to be a number or string, not a ',class(seed)])
end

% and now the CGR iteration
M(:,1)=CGRo+(M(:,1)-CGRo).*0.5;
for i=2:m
    M(:,i)=M(:,i-1)+(M(:,i)-M(:,i-1)).*0.5;
end


% ----------- NESTED FUNCTION ----------------
% ----------- kernel functions ---------------
function k=usm_kernel(u,L,S)

%USM_KERNEL determines density distribution from USM coordinates
%Syntax: k=usm_kernel(u,L,S)
%Description: u is a set of UNIDIRECTIONAL usm coordinates (~ the same as
%             CGR coordinates), L is the memory length and S is the
%             smoothing

[D,n]=size(u);
% Set centroids
ui=[1/2^(L+1):1/2^L:1];
% Set kernel size
un=num2str(length(ui));uk='';for i=1:D;uk=[uk,',',un];end;uk=str2num(uk(2:end));
if length(uk)==1
    k=zeros(1,uk); %remember zeros(x) is the same as zeros(x,x) when x is a scalar (that was silly Mathworks!)
else
    k=zeros(uk);
end
nn=prod(uk); %number elements in k
Ind=[1:nn];co=cell(1,D);[co{:}]=ind2sub(uk,Ind);co=cell2mat(co');% generate all coordinates
%SS=sum(S.^[0:L]);
ff=find(u~=1);%in some very rare cases u may have value 1 so that needs not to be pushed past the last quadrant
for i=0:L %for each memory length
    H=(((2^D)*S)^i)/sum(S.^[0:L]); %calculate height for (i+1) order, = ith memory length
    %H=(((2^D)/S)^i)/sum(S.^(-[0:L])); %<-- replace S with 1/S after manuscript is accepted, it makes more sense
    %Apply it to all the relevant quandrants
    qmin=floor(u*(2^i));qmin(ff)=qmin(ff)*2^(L-i)+1; %see note about ff earlier
    qmax=qmin+2^(L-i)-1;
    %qmax=ceil(u*(2^i));qmax=qmax*2^(L-i);
    for j=1:n %for each unit of the sequence
        f=find(prod((co>=repmat(qmin(:,j),1,nn)).*(co<=repmat(qmax(:,j),1,nn)),1));
        k(f)=k(f)+H;
    end
end