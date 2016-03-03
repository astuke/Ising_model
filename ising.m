function ising()
    
    tic;
    
    L=8; %length of the lattice
    p=linspace(0.0,1.0,20); %probability that the coupling constant J is -1 for a nearest neighbour pair,
                            % see function J=RandomCoupling(p) below
    
    NTimes=100;
    
    M=zeros(NTimes,length(p)); % Magnetization
    Ms=zeros(NTimes,length(p)); % Staggered Magnetization

    %Calculate the expectation value of M and Ms as a function of p
    parfor k=1:length(p)
        
        [M(:,k),Ms(:,k)]=AnnealMultipleTimes(p(k),L,NTimes)
        
    end
    
    m=mean(M); %expectation value of M
    ms=mean(Ms); %expectation value of Ms
    ErrM=std(M)/sqrt(NTimes); %error estimate of expectation value (std: standard deviation)
    ErrMs=std(Ms)/sqrt(NTimes);
    
    figure;
    hold on;
    errorbar(p,m,ErrM,'b');
    errorbar(p,ms,ErrMs,'r');
    xlabel('p');
    legend('m','m_s');
    title(['L=',num2str(L)]);
    hold off;
    
    figure;
    hold on;
    plot(p,ErrM,'b*');
    plot(p,ErrMs,'r*');
    hold off;
    
    toc;
    [couplings,J]=CreateCouplingTables(L,p);
    [M,Ms,spins]=DoAnnealing(p,L);
    figure;
    PlotCouplingTables(couplings,J,L);
    PlotSpins(spins,L);
        
end

%This function calculates the ground state N times for different realizations of the coupling constants, 
%obtaining the value of M and Ms for each ground state
function [M,Ms]=AnnealMultipleTimes(p,L,NTimes) 
    M=zeros(NTimes,1);
    Ms=zeros(NTimes,1);

    for k=1:NTimes
        [M(k),Ms(k)]=DoAnnealing(p,L);
    end

end

%This is the simulated annealing procedure to calculate the ground state of the system.  
function [M,Ms,spins]=DoAnnealing(p,L)
    
    [couplings,J]=CreateCouplingTables(L,p);
    
    NUpdates=1000;
    
    T=linspace(2.0,0.2,NUpdates); %temperature
    beta=1./T; %inverse temperature
    
    spins=randi([0,1],1,L*L)*2-1;
    
    %Start running the Monte Carlo algorithm (UpdateSweep) at a high temperature and then cool down
    for k=1:NUpdates
        spins=UpdateSweep(beta(k),spins,couplings,J);
    end
    
    M=magnetization(spins,L);
    Ms=StaggeredMagnetization(spins,L);
    
end

function Ms=StaggeredMagnetization(spins,L)

    Ms=0;

    for k=1:length(spins)
        [x,y]=ExpandIndex(k,L);
        sign=(-1)^(x+y);
        Ms=Ms+sign*spins(k);
    end
    
    Ms=abs(Ms)/(L*L);

end

function M=magnetization(spins,L)

    M=abs(sum(spins))/(L*L);
    
end

%Monte Carlo algorithm: Spin flip is accepted with certain probability
function spins=UpdateSweep(beta,spins,couplings,J)

    N=length(spins);
    
    for k=1:length(spins)
        site=randi([1,N]);
        %propose to flip the chosen spin
        DeltaE=2*spins(site)*sum(spins(couplings(site,:)).*J(site,:)); % energy change of system by flipping spin
        
        %accept spin flip with probabilty
        if DeltaE<0 || rand()<exp(-beta*DeltaE)
            spins(site)=-spins(site); %update spin
        end
            
    end
    
end

%Plot the spins given by the table spins. The spins are given as a single
%vector whose length is L*L. spins(I) then corresponds to the site at
%coordinates [x,y]=ExpandIndex(I,L).
function PlotSpins(spins,L)

    hold on;
    
    for k=1:length(spins)
       
        [x,y]=ExpandIndex(k,L);
        
        if spins(k)>0
            plot(x,y,'r*');
        else
            plot(x,y,'b*');
        end
        
    end
    
    hold off;

end

%This function plots the coupling tables created by CreateCouplingTables.
function PlotCouplingTables(couplings,J,L)

    hold on;
    for k=1:size(couplings,1)
        
        [x,y]=ExpandIndex(k,L);
        
        %This plots the couplings to "down" and "right" directions.
        for l=1:2
            
            [xn,yn]=ExpandIndex(couplings(k,l),L);
            
            %no over boundary hoppings plotted
            if abs(xn-x)<2 && abs(yn-y)<2
                if J(k,l)>0
                    plot([x,xn],[y,yn],'r-');
                else
                    plot([x,xn],[y,yn],'b-');
                end
            end
            
        end
         
        %This plots also the "up" and "left" directions. 
        for l=3:4
           
           [xn,yn]=ExpandIndex(couplings(k,l),L);
           
           %no over boundary hoppings plotted
           if abs(xn-x)<2 && abs(yn-y)<2
               if J(k,l)>0
                   plot([x,xn],[y,yn],'r--');
               else
                   plot([x,xn],[y,yn],'b--');
               end
           end
        end    
        
        
    end
    hold off;

end

%Creates two tables, couplings and J. The couplings table is a (L*L)x4
%matrix, which tells which sites are neighbours of a given site. Each site
%has a unique index that can be calculated from the corresponding
%coordinates using the function FlattenIndex. The neighbours of site number
%I are then given by couplings(I,:), so that if e.g.
%couplings(I,:)=[1,2,3,4], then the sites 1,2,3 and 4 would be neighbours
%of I. The corresponding coupling constants are then given by the J matrix
%as e.g. J(I,:)=[1,1,-1,1].
function [couplings,J]=CreateCouplingTables(L,p)
    
    couplings=zeros(L*L,4);
    J=zeros(L*L,4);
    
    for x=1:L
        for y=1:L
            
            I=FlattenIndex(x,y,L);
            
            xright=WrapAdd(x,1,L);
            ydown=WrapAdd(y,-1,L);
            
            IRight=FlattenIndex(xright,y,L);
            IDown=FlattenIndex(x,ydown,L);
            
            JRight=RandomCoupling(p);
            JDown=RandomCoupling(p);
            
            %couplings right, down, left, up
            couplings(I,1)=IRight;
            couplings(I,2)=IDown;
            couplings(IRight,3)=I;
            couplings(IDown,4)=I;
            
            J(I,1)=JRight;
            J(I,2)=JDown;
            J(IRight,3)=JRight;
            J(IDown,4)=JDown;
            
        end
    end
    
end

%Calculate the index of the site corresponding to the 2d coordinates [x,y].
%L is the size of the lattice. (This is the "inverse function" of
%ExpandIndex.)
function I=FlattenIndex(x,y,L)
    
    I=(x-1)+(y-1)*L+1;
    
end

%Calculate the coordinates corresponding to site number I. L is the size of
%the LxL lattice. (This is the "inverse function" of FlattenIndex.)
function [x,y]=ExpandIndex(I,L)

    x=mod(I-1,L)+1;
    y=fix((I-1)/L)+1;
    
end

%Move the coordinate x by a distance dx taking into account periodic
%boundary conditions.  L is the size of the LxL lattice.
function xnew=WrapAdd(x,dx,L)

    xnew=mod(x-1+dx,L)+1;

end

%Randomly choose the coupling J to be +-1 based on the probability p. 
function J=RandomCoupling(p)
    
    if rand()<p
        J=-1; %negative coupling: neighbouring soins antiparallel
    else
        J=1; %positive coupling: neighbouring spins parallel. For example, if p=0, it is J=1 for all nearest neighbour pairs,
             %so that all spins on the lattice are aligned (all up or all down). This is the completely magnetized ground state.
    end

end
