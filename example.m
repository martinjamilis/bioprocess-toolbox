clear, clc, close all

a=BioProcess(3,2)
a.RateArray{1,1}=kineticModel('proportional',1);
a.RateArray{2,1}=kineticModel('proportional',1);
a.RateArray{1,2}=kineticModel('monod',0.5,0.5);
a.RateArray{1,3}=kineticModel('inhibit',1,10);
a.RateArray{2,2}=kineticModel('haldane',0.7,0.6,5);

a.YieldMatrix = [1 0;-0.2 -0.3;0 1];
xi_in=[0;10;0];


a.stateSpace(1,[1;2;3],[0;10;0],0.1)