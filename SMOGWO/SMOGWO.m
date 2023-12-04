clc;
clear;
close all;

TestProblem='ZDT1';
nVar=10;

fobj = cec09(TestProblem);

xrange = xboundary(TestProblem, nVar);

% Lower bound and upper bound
lb=xrange(:,1)';
ub=xrange(:,2)';

VarSize=[1 nVar];

GreyWolves_num=100;
MaxIt=20;  % Maximum Number of Iterations
Archive_size=100;   % Repository Size
Ks=50;     % 邻域圈的大小

runtimes=3;
GD=zeros(runtimes,1);
IGD=zeros(runtimes,1);
Spread=zeros(runtimes,1);
TPF=xlsread('Case3.xls',TestProblem);

alpha=0.1;  % Grid Inflation Parameter
nGrid=10;   % Number of Grids per each Dimension
beta=4; %=4;    % Leader Selection Pressure Parameter
gamma=2;    % Extra (to be deleted) Repository Member Selection Pressure

% Initialization

GreyWolves=CreateEmptyParticle(GreyWolves_num);
X=CreateEmptyParticle(3*GreyWolves_num);    

for r=1:runtimes

    
for i=1:GreyWolves_num
    GreyWolves(i).Velocity=0;
    GreyWolves(i).Position=zeros(1,nVar);
    for j=1:nVar
        GreyWolves(i).Position(1,j)=unifrnd(lb(j),ub(j),1);
    end
    GreyWolves(i).Cost=fobj(GreyWolves(i).Position')';
    %%%%%%%%%%%反向学习策略
    position1=abs(rand(1,nVar).*(lb+ub)-GreyWolves(i).Position);
    position1=min(max(position1,lb),ub);
    position1_value=fobj(position1')';
    if Dominates(position1_value,GreyWolves(i).Cost)
        GreyWolves(i).Position=position1;
        GreyWolves(i).Cost=position1_value;
    end
    %%%%%%%%%%%%
    GreyWolves(i).Best.Position=GreyWolves(i).Position;
    GreyWolves(i).Best.Cost=GreyWolves(i).Cost;
end

GreyWolves=DetermineDomination(GreyWolves);

Archive=GetNonDominatedParticles(GreyWolves);

Archive_costs=GetCosts(Archive);
G=CreateHypercubes(Archive_costs,nGrid,alpha);

for i=1:numel(Archive)
    [Archive(i).GridIndex, Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
end

% MOGWO main loop

for it=1:MaxIt
    a=2-2*(it/MaxIt);
    for i=1:GreyWolves_num
        
%         clear rep2
%         clear rep3
        
        % Choose the alpha, beta, and delta grey wolves
        Delta=SelectLeader(Archive,beta);
        Beta=SelectLeader(Archive,beta);
        Alpha=SelectLeader(Archive,beta);
        
        % If there are less than three solutions in the least crowded
        % hypercube, the second least crowded hypercube is also found
        % to choose other leaders from.
        if size(Archive,1)>1
            counter=0;
            for newi=1:size(Archive,1)
                if sum(Delta.Position~=Archive(newi).Position)~=0
                    counter=counter+1;
                    rep2(counter,1)=Archive(newi);
                end
            end
            Beta=SelectLeader(rep2,beta);
        end
        
        % This scenario is the same if the second least crowded hypercube
        % has one solution, so the delta leader should be chosen from the
        % third least crowded hypercube.
        if size(Archive,1)>2
            counter=0;
            for newi=1:size(rep2,1)
                if sum(Beta.Position~=rep2(newi).Position)~=0
                    counter=counter+1;
                    rep3(counter,1)=rep2(newi);
                end
            end
            Alpha=SelectLeader(rep3,beta);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%选择领域
        Q11=randperm(GreyWolves_num,Ks);
       GreyWolves11=GreyWolves(Q11);
       GreyWolves11=DetermineDomination(GreyWolves11);
       
       Archive11=GetNonDominatedParticles(GreyWolves11);
       
       a1=length(Archive11);
       k11=randperm(a1,1);
       
       moution11= Archive11(k11).Position; % 邻域内非支配解
        
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      改进2
      %%% 头狼变异策略
        moution1=Delta.Position;
        r1=randperm(nVar,1); 
        kk=randperm(nVar,1);
        while(r1==i)
            r1=randperm(nVar,1);
        end
        moution2=GreyWolves(r1).Position;
        moution1(kk)=moution11(kk)+(moution11(kk)-moution2(kk))*(rand-0.5)*2+(moution1(kk)-moution11(kk))*rand*1.5;
        moution1=min(max(moution1,lb),ub);
        moution1_value=fobj(moution1')';
        if Dominates(moution1_value,Delta.Cost)
            Delta.Position=moution1;
            Delta.Cost=moution1_value;
        end
        
        moution3=Beta.Position;
        r2=randperm(nVar,1); 
        kk=randperm(nVar,1);
        while(r2==i)
            r2=randperm(nVar,1);
        end
        moution4=GreyWolves(r2).Position;
        moution3(kk)=moution11(kk)+(moution11(kk)-moution4(kk))*(rand-0.5)*2+(moution3(kk)-moution11(kk))*rand*1.5;
        moution3=min(max(moution3,lb),ub);
        moution3_value1=fobj(moution3')';
        if Dominates(moution3_value1,Beta.Cost)
            Beta.Position=moution3;
            Beta.Cost=moution3_value1;
        end
        
        moution5=Alpha.Position;
        r3=randperm(nVar,1); 
        kk=randperm(nVar,1);
        while(r3==i)
            r3=randperm(nVar,1);
        end
        moution6=GreyWolves(r3).Position;
        moution5(kk)=moution11(kk)+(moution11(kk)-moution6(kk))*(rand-0.5)*2+(moution5(kk)-moution11(kk))*rand*1.5;
        moution5=min(max(moution5,lb),ub);
        moution5_value2=fobj(moution5')';
        if Dominates(moution5_value2,Alpha.Cost)
            Alpha.Position=moution5;
            Alpha.Cost=moution5_value2;
        end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        % Eq.(3.4) in the paper
        c=2.*rand(1, nVar);
        % Eq.(3.1) in the paper
        D=abs(c.*Delta.Position-GreyWolves(i).Position);
        % Eq.(3.3) in the paper
        A=2.*a.*rand(1, nVar)-a;
        % Eq.(3.8) in the paper
        X1=Delta.Position-A.*abs(D);
        
        
        % Eq.(3.4) in the paper
        c=2.*rand(1, nVar);
        % Eq.(3.1) in the paper
        D=abs(c.*Beta.Position-GreyWolves(i).Position);
        % Eq.(3.3) in the paper
        A=2.*a.*rand()-a;
        % Eq.(3.9) in the paper
        X2=Beta.Position-A.*abs(D);
        
        
        % Eq.(3.4) in the paper
        c=2.*rand(1, nVar);
        % Eq.(3.1) in the paper
        D=abs(c.*Alpha.Position-GreyWolves(i).Position);
        % Eq.(3.3) in the paper
        A=2.*a.*rand()-a;
        % Eq.(3.10) in the paper
        X3=Alpha.Position-A.*abs(D);
        
        % Eq.(3.11) in the paper
        GreyWolves(i).Position=(X1+X2+X3)./3;
        
         % Boundary checking
        X111=min(max(X1,lb),ub);
        X222=min(max(X2,lb),ub);
        X333=min(max(X3,lb),ub);
        
        X(i).Position=X111;
        X(GreyWolves_num+i).Position=X222;
        X(2*GreyWolves_num+i).Position=X333;
        GreyWolves(i).Position=min(max(GreyWolves(i).Position,lb),ub);
        
        
        X(i).Cost=fobj(X111')';
        X(GreyWolves_num+i).Cost=fobj(X222')';
        X(2*GreyWolves_num+i).Cost=fobj(X333')';
        GreyWolves(i).Cost=fobj(GreyWolves(i).Position')';
 
    end
    
     GreyWolves=DetermineDomination(GreyWolves);
     X=DetermineDomination(X);
     non_dominated_wolves=GetNonDominatedParticles(GreyWolves);
     non_dominated_X=GetNonDominatedParticles(X);
    
    temp=[Archive
        non_dominated_wolves
        non_dominated_X];
    
    Archive=DetermineDomination(temp);
    Archive=GetNonDominatedParticles(Archive);
    
    for i=1:numel(Archive)
        [Archive(i).GridIndex, Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
    end
    
    if numel(Archive)>Archive_size
        EXTRA=numel(Archive)-Archive_size;
        Archive=DeleteFromRep(Archive,EXTRA,gamma);
        
        Archive_costs=GetCosts(Archive);
        G=CreateHypercubes(Archive_costs,nGrid,alpha);
        
    end
    
%     disp(['In iteration ' num2str(it) ': Number of solutions in the archive = ' num2str(numel(Archive))]);
%     save results
    
    % Results
    
    costs=GetCosts(GreyWolves);
    Archive_costs=GetCosts(Archive);
    
%     if drawing_flag==1
%         hold off
%         plot(costs(1,:),costs(2,:),'k.');
%         hold on
%         plot(Archive_costs(1,:),Archive_costs(2,:),'rd');
%         legend('Grey wolves','Non-dominated solutions');
%         drawnow
%     end
    
end
plot(Archive_costs(1,:),Archive_costs(2,:),'rO');
GD(r)=GD_matlab(Archive_costs',TPF);
IGD(r)=IGD_matlab(Archive_costs',TPF);
Spread(r)=Spread_matlab(Archive_costs',TPF);
end
mean(GD)
mean(IGD)
mean(Spread)
aa=Archive_costs';
