%Author Pouya Tavakoli
%Article: Generating synthetic ground motions reaching target spectrum with
%the optimization approach (2023)
clear;clc;close all;
% in this file the optimization algorithm is implemented
strfitnessfct = 'ObjFF'; %the objective funtion
N=200; % number of optimization variables
%The numbers 1 to 100 of the optimization vatiables are used for determining the
%amplitude of the earthquake, and the numbers 101 to 200 are used to
%determinig the phase angles

%% determining the boundaries of the optimization variables
%Boundaries of the amplitudes
for i=1:100;
R(1:2,i) = [0 0.07]';    % Bounds on decision variables 
end
%according this code, the optimization algorithm has to select the values
%in the range of zero to 0.07 in the scale of g (gravity of earth). To put
%it in another words, if for a given optimization value, the 0.07 is
%selected from the optimization algorithm, in the corresponding frequency
%the amplitude of the sinusoidal function will be 0.07g. 

%Boundaries of the phase angles
for i=101:200;
    R(1:2,i) = [0 2*pi]';
end
%% CMAES Algorithm 
%for more information readers can refer to Hansen N. The CMA evolution strategy: A tutorial. arXiv preprint arXiv:
%1604.00772. 2016.
lb = R(1,:)';ub = R(2,:)';xmean=rand(N,1).*(ub - lb) + lb;
sigma = 0.3;stopfitness = 0.01;     lambda =20       ;stopeval =1500* lambda;
% in the stopeval, the number of 1500 is about the number of the iterations.
mu = lambda/2;weights = log(mu+1/2)-log(1:mu)';mu = floor(mu);     
weights = weights/sum(weights);  mueff=sum(weights)^2/sum(weights.^2);
cc = (4+mueff/N) / (N+4 + 2*mueff/N);cs = (mueff+2) / (N+mueff+5);
c1 = 2 / ((N+1.3)^2+mueff);cmu = min(1-c1, 2 * (mueff-2+1/mueff) / ((N+2)^2+mueff));
damps = 1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs; 
pc = zeros(N,1); ps = zeros(N,1);B = eye(N,N);Dd = ones(N,1);C = B * diag(Dd.^2) * B';
invsqrtC = B * diag(Dd.^-1) * B';eigeneval = 0;chiN=N^0.5*(1-1/(4*N)+1/(21*N^2)); 
counteval = 0;
it=1;

while counteval < stopeval
    for k=1:lambda
        arx(:,k) = xmean + sigma * B * (Dd .* randn(N,1)); 
        arx(:,k) = max(arx(:,k),lb);arx(:,k) = min(arx(:,k),ub);
        arfitness(k) = feval(strfitnessfct, arx(:,k)');  counteval = counteval+1;
    end
  
    [arfitness, arindex] = sort(arfitness);arfitness(1);
    Best(it)=arfitness(1);
    Best2(it)=min(Best);
%the objective funciton iteration curve is plotted here in each iteration.
      plot(Best2,'b-') ; title(['Best = ', num2str(Best2(it))]) 
     xlabel('Iteration')
     ylabel('Objective function') ; drawnow
      it=it+1;
    xold = xmean;xmean = arx(:,arindex(1:mu))*weights;   
    ps = (1-cs)*ps + sqrt(cs*(2-cs)*mueff) * invsqrtC * (xmean-xold) / sigma; 
    hsig = norm(ps)/sqrt(1-(1-cs)^(2*counteval/lambda))/chiN < 1.4 + 2/(N+1);
    pc = (1-cc)*pc + hsig * sqrt(cc*(2-cc)*mueff) * (xmean-xold) / sigma;
    artmp = (1/sigma) * (arx(:,arindex(1:mu))-repmat(xold,1,mu));
    C = (1-c1-cmu) * C + c1 * (pc*pc'+ (1-hsig) * cc*(2-cc) * C) ...
         + cmu * artmp * diag(weights) * artmp'; 
    sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1));
    if counteval - eigeneval > lambda/(c1+cmu)/N/10
        eigeneval = counteval;C = triu(C) + triu(C,1)';[B,Dd] = eig(C);
        Dd = sqrt(diag(Dd));invsqrtC = B * diag(Dd.^-1) * B';
    end
    if arfitness(1) <= stopfitness || max(Dd) > 1e7 * min(Dd); break; end
    
end
x = arx(:, arindex(1)),fval=arfitness(1)

save xmin.txt x -ascii % the optmimum values of the 200 variables are saved in the Xmin file 
% this file is imported to the Plott file so as to extract the detail
% information of the earthquake.
save Best_sandiegoAspectra_CMAES.txt Best2 -ascii %the values of the objective function is saved in each iteration