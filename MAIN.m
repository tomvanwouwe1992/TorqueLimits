% This is the header file that calls the TorqueLimitClassifier function

% Create some dummy data to test
size = 100;   % number of samples
q = 0.3*rand(size,5) - 0.15*rand(size,5);   % Select some realistic value for q and qdot
q_dot = 1*rand(size,5) - 0.5*rand(size,5);
tau = [50*rand(size,3) - 25*rand(size,3), 100*rand(size,2) - 50*rand(size,2)];  % Select some realistic value for the torques
% These dummy values seem to create a nice mix of feasible and unfeasible
% torques

% Call the classifier, it returns a vector which indicates whether the
% chosen sample is feasible, activations that are needed and the values of
% the reserveactuators.
[Torque_Manageable,Activations,ReserveActuators] = TorqueLimitClassifier(q,q_dot,tau);

% Note: In case of infeasible torque (indicated by 0) the activations
% should max out, i.e. reach 1 (at least for some muscles) and the reserve
% actuator torques are non-trivial.