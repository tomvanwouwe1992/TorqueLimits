size = 100;
q = 0.3*rand(size,5) - 0.15*rand(size,5);
q_dot = 1*rand(size,5) - 0.5*rand(size,5);
tau = [50*rand(size,3) - 25*rand(size,3), 100*rand(size,2) - 50*rand(size,2)]; 
[Torque_Manageable,Activations,ReserveActuators] = TorqueLimitClassifier(q,q_dot,tau);