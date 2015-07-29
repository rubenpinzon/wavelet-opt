%Script that demostrates the procedure to optimize a wavelet function using
%GA and the lifting scheme following the procedure introduced in [1,2].
%
%[1]R D Pinzon-Morales and Yutaka Hirata. 2012. Customization of
%wavelet function for pupil fluctuation analysis to evaluate levels of 
%sleepiness. In Proceedings of the 11th international conference on Teleco
%mmunications and Informatics, Proceedings of the 11th international 
%conference on Signal Processing (SITE'12),115-120. 
%[2]R D Pinzon Morales et al. Novel signal-dependent filter bank method for
%identification of multiple basal ganglia nuclei in Parkinsonian patients
%2011 J. Neural Eng. 8 036026

%Step 1: Extract from the database a training set for the algorithm. 
%        If the database is labeled, then the training set should contain
%        samples from each class. For the example we will create 100
%        signals with 1000 samples from gaussians with muA = 2,
%        sigmaA = 0.1, and for classB, muB = 1, sigma = 0.9

classA  = normrnd(2,0.1,1000,100);
classB  = normrnd(1,0.9,1000,100);

y_train = [classA, classB];
labels  = [ones(1,100), 2*ones(1,100)];

%Step 2: Setup the options for the GA procedure. The GA is implemented with
%        the library GAOT from Christopher R Houck: chouck@eos.ncsu.edu

N           = 6;   %Order of the wavelet filters which is the same as the 
                   %order of the FIR filters: Update and Predictor in the
                   %implementation of the WT using Lifting Schemes
level       = 5;   %Wavelet decomposition level: Application depended as N

costFnc     = 'MyCostFunc'; %name of the m-file containing the cost functi
                            %on to be maximized, for example, class
                            %separability

options{1}  = N;                   %filter order  
options{2}  = level;               %decomposition level
options{3}  = y_train;             %training data
options{4}  = labels;              %labels

numPop      = 20;                  %population scale GA 
termOps     = 10;                  %Generations of the GA
nu_ga       = 10;                  %How many repetitions of the whole GA to
                                   %look up for local solutions
bounds      = repmat([0.25 -0.25],N-2,1); %P y U bounds N-2
vals        = 0;

%Step 3: Main cicle of the optimization

for gen=1:nu_ga
    % Initial population of the GA
    initPop = initializega(numPop,bounds,costFnc,options);
    %Evolve the initial population following GA rules. Parameters of the Ga
    %can be consuted in the Manual of the GAOT Toolbox. By default:
    
    [p, endPop, bestSols,trace]=ga(bounds,'objEvalFMI',options,initPop,...
        [1e-6 1 1],'maxGenTerm',termOps,'normGeomSelect',0.08,...
        'arithXover',[2 0],'nonUnifMutation',[2 1 3]);
    
    %reconstruct the wavelet operators P and U from the GA results
    P=0.5-sum(p(1:floor(N/2)-1))*[1 1];      %Po=P1=sum(p) Filter normalization
    U=0.25-sum(p(floor(N/2):end-1))*[1 1];
    for i=floor(N/2)-1:-1:1                  %Symmetrical linear phase constraint 
        P=[p(i) P p(i)];
        U=[p(i+floor(N/2)-1) U p(i+floor(N/2)-1)];
    end
    
    %Track the resutls 
    Uf(gen,:)=U;
    Pf(gen,:)=P;
end

%Step 4(optional):Visualize the wavelet funtions evolved using matlab
%wavelet toolbox
for i=1:nu_ga
    P           = Pf(i,:);
    U           = Uf(i,:);
    LS          = cell(3,3);
    LS{1,1}     ='d';
    LS{1,2}     =(-1)*P; 
    LS{1,3}     =[floor(N/2)];
    LS{2,1}     ='p';
    LS{2,2}     = U; 
    LS{2,3}     =[floor(N/2)];
    LS{3,1}     =[1];
    LS{3,2}     =[1];
    [LoD,HiD]   =ls2filters(LS,'a_num');   
    [vx,phi(:,i),psi(:,i)]  =   phivals(LoD',3);
    h0          =LoD;
    h1          =HiD;
    N           = 512;
    W           = 2/N*(-N/2:N/2-1);
    H0          = fftshift(fft(h0,N));
    H1          = fftshift(fft(h1,N));
    hold on;
    plot(W, abs(H0)/max(abs(H0)), '--','DisplayName','Low pass filter H_0'),hold on
    plot(W, abs(H1),'DisplayName','High pass filter H_1')
    xlabel('Angular frequency (normalized by \pi)')
    ylabel('Frequency response magnitude')
    filtros(i,:)=abs(H0)/max(abs(H0));    
end
errorbar(0:80,mean(phi'),std(phi')), hold on
errorbar(0:80,mean(psi'),std(psi')), hold on
figure(2)
plot(mean(phi')),plot(mean(psi'))
errorbar(mean(filtros),std(filtros)), hold on
plot(mean(filtros))

