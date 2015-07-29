function [sol, val] = MaCostFunc(sol,options,ctr)
%This function contains the desired metric to optimize the wavlet function
%the in/out parameters must remain unchanged, otherwise Initialize.m has to
%be changed accordingly.
%
%Inputs: sol (1 x N-2) : Current values of the filters U and P
%        options (cell): contains options passed from the calling function 
%                        passing throuhg the initializega.m
%        options{1}    : N, filters U and P order
%        options{2}    : wavelet decomposition level
%        options{3}    : X, trainig database
%        options{4}    : label, class labels of the training database
%        ctr           : counter to track the progress   
%Output: sol
%        val           : value of the cost function after processing the
%                        current individual. The GA maximize this value.

tic
%
%Ruben Dario Pinzon Morales

N       = options{1};               %filter's order
level   = options{2};               %label vector 
X       = options{3};               %training data (samples X signals)
label   = options{4};               %Labels (1 x signals) 
p       = sol(1:floor(N/2)-1);      %predictor's coefficients
u       = sol(floor(N/2):end-1);    %update coefficients

P=[0.5-sum(p(1:end-1)) 0.5-sum(p(1:end-1))];%Po=P1 Filter normalization linear constraint
U=[0.25-sum(u(1:end-1)) 0.25-sum(u(1:end-1))];%Uo=U1 Filter normalization linear constraint

for i=floor(N/2)-1:-1:1                 %Symmetrical linear phase constraint 
    P=[p(i) P p(i)]; %#ok<AGROW>
    U=[u(i) U u(i)]; %#ok<AGROW>
end

%% wavelet decomposition and cost function computation (class separatibilty as example)
numsigs=size(X,1);
numNodes=2^(level+1)-1;

%Descompose the trainig signals using the current individuals U and P and
%the wavlet packet
W=zeros(numNodes,numsigs);       %for wavelet energy decomposition 
for s=1:size(X,1)      
    W(:,s)=wplift(X(s,:),U,P,level); 
end
%Calculate Class Separability Index of the wavelet packet coefficients
for i=1:max(label)
    Ni=sum(label==i);    
    for k=1:numNodes %C
        y_c=W(k,label==i);
        m_c=mean(y_c);
        dii(k)=(1/Ni*sum((y_c-m_c).^2)^(1/2));
        mt(i,k)=mean(y_c);
    end
    Dii(i)=sum(dii.^2)^(1/2);
end
Mii=mean(mt,2);
comb=max(label)-1;
i=0;
Rij=[];
while (comb~=0)
    i=i+1;
    ind=0;
    for j=i+1:max(label)
        ind=ind+1;
        rij(ind)=(Dii(i)+Dii(j))./abs(Mii(i)-Mii(j));
    end
    Rij=[Rij rij];
    comb=comb-1;
    clear rij
end
%=========================IMPORTANT===========================
val=min(Rij);%<<<<<<<<This is the value to Maximize. Change according to 
             %the specific objectives of your application 
             
fprintf('Gen -> %d-th, Took: %3.1f s, Cost value %4.4f\n',ctr,toc,val)
