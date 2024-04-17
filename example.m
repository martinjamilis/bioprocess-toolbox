clear, clc, close all

% s->x+p
a=BioProcess(3,1)
a.Names={'x', 's', 'p'};
a.RateArray{1,1}=kineticModel('proportional',1);
a.RateArray{1,2}=kineticModel('monod',0.5,0.5);
a.YieldMatrix = [1;-1/0.5;0.5];
xi_in=[0;10;0];

% s=linspace(0,10,100);
% for idx=1:100
%     mu(idx)=a.rates([1;s(idx)]);
% end
% figure(1); plot(s,mu)

[t,xi]=a.euler(0.001,[0 48],[0.1;20;0],[0;20;0],0.1);
figure(2), plot(t,xi)

%% test
clear, clc, close all

% x->s
a=BioProcess(3,2)

a.Names={'x', 's', 'p'};

% a.RateArray{1,1}=kineticModel('proportional',1);
a.RateArray{1,2}=kineticModel('monod',0.5,0.5);
% a.RateArray{2,1}=kineticModel('proportional',1);
% a.RateArray{1,3}=kineticModel('inhibit',1,10);
% a.RateArray{2,2}=kineticModel('haldane',0.7,0.6,5);

% a.RateArray{2,2}=5;

a.YieldMatrix = [1 0;-0.2 -0.3;0 1];
xi_in=[0;10;0];


a.stateSpace(1,[1;2;3],[0;10;0],0.1)

[t,xi]=a.euler(0.01,[0 48],[0.1;1;0],[0;0;0],0.25);
figure, plot(t,xi)