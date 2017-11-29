function [Torque_Manageable,act,eT] = TorqueLimitClassifier(q, qdot, tau)
% This function calculates fot he lower limb whether it can produce the
% provided torque under the specified kinematic state of the model.

import org.opensim.modeling.*
MAIN_path = pwd;

load ActiveFVParameters
load PassiveFLParameters
load Faparam

N = size(q,1);


% Load the OpenSimModel
Model_OS = Model(fullfile(MAIN_path,'gait2354_simbody.osim'));
DOF = {'hip_flexion_l' 'hip_adduction_l' 'hip_rotation_l' 'knee_angle_l' 'ankle_angle_l'};
nDof = length(DOF);
[MusclesOI,CoordinatesOI] = getMuscleStructure(Model_OS,DOF);
M = length(MusclesOI);
act = ones(N,M);
FMltilde = ones(N,M);
FMvtilde = ones(N,M);
Fpe = ones(N,M);
cos_alpha = ones(N,M);
% Translate Coordinate information into DOF_number
state = [q, qdot];
[LMT,VMT,dM] = get_LMT_vMT_dM(Model_OS,MusclesOI,CoordinatesOI,state,DOF);

for m = 1:M
    [~, ~, FMltilde(:,m), FMvtilde(:,m), Fpe(:,m), cos_alpha(:,m)] = HillModel_RigidTendon(act(:,m),LMT(:,m),VMT(:,m),MusclesOI(m),ActiveFVParameters,PassiveFLParameters,Faparam);
end


% Equations TBC
FMo = ones(size(act,1),1)*[MusclesOI(:).maxIsoForce];
Fpas = FMo.*Fpe.*cos_alpha;
Fact = FMo.*FMltilde.*FMvtilde.*cos_alpha;


% Add optimal force for reserve torques
Topt = 1;
% Topt = 150;
Fpas = [Fpas zeros(N,nDof)];
Fact = [Fact Topt*ones(N,nDof)];

ID_data = tau;

x0 = repmat(0.01*ones(M+nDof,1),N,1);


% The bounds
options.lb = repmat([zeros(M,1); -1500*ones(nDof,1)],N,1);
options.ub = repmat([ones(M,1); 1500*ones(nDof,1)],N,1);
options.cl = repmat(zeros(nDof,1),N,1);
options.cu = repmat(zeros(nDof,1),N,1);

I = N*(M+nDof);
% Set up the auxiliary data.
options.auxdata = { M N nDof reshape(Fact', I, 1)  reshape(Fpas', I, 1) ...
    ID_data dM};


options.ipopt.mu_strategy      = 'adaptive';
options.ipopt.max_iter         = 1500;
options.ipopt.tol              = 1e-5;
options.ipopt.hessian_approximation = 'limited-memory';

% The callback functions.
funcs.objective         = @objective_SO;
funcs.constraints       = @constraints_SO;
funcs.gradient          = @gradient_SO;
funcs.jacobian          = @jacobian_SO;
funcs.jacobianstructure = @jacobianstructure_SO;

[x,~] = ipopt_auxdata(x0,funcs,options);

x_opt = reshape(x, M+nDof, N)';

act = x_opt(:,1:M);
eT = x_opt(:, M+1:M+nDof)*Topt;

ReserveActuatorSum = sum(abs(eT),2);
for i = 1:size(ReserveActuatorSum,1)
    if ReserveActuatorSum(i) < 1
        ReserveActuatorSum(i) = 1;
    else 
        ReserveActuatorSum(i) = 0;
    end
end


Torque_Manageable = ReserveActuatorSum;
end
% ------------------------------------------------------------------
function f = objective_SO (x, auxdata)
%f     = 0.5 * sum(x.^2);
[M, N, nDof, Fmax, Fpas, ID_data, MomentArm] = deal(auxdata{:});
x_opt = reshape(x, M+nDof, N)';
f     = 0.5 * (sum(sum(x_opt(:,1:M).^2)) + sum(sum(x_opt(:,M+1:end).^2)));
end

% ------------------------------------------------------------------
function c = constraints_SO (x, auxdata)
[M, N, nDof, Fmax, Fpas, ID_data, MomentArm] = deal(auxdata{:});

F = Fmax .* x + Fpas;

c = zeros(nDof*N,1);

for k = 1:nDof
    F_matrix = reshape(F, M+nDof, N)';
    MomentArm_matrix = reshape(MomentArm(:,k), M+nDof, N)';
    c(k:nDof:end) = sum(F_matrix.*MomentArm_matrix, 2) - ID_data(:,k);
end
end
% ------------------------------------------------------------------
function g = gradient_SO (x, auxdata)
%g = x;
[M, N, nDof, Fmax, Fpas, ID_data, MomentArm] = deal(auxdata{:});
x_opt = reshape(x, M+nDof, N)';
gtemp = [x_opt(:,1:M) x_opt(:,M+1:end)];
I = N*(M+nDof);
g = reshape(gtemp',I,1);
end
% ------------------------------------------------------------------
function J = jacobianstructure_SO (auxdata)
[M, N, nDof, Fmax, Fpas, ID_data, MomentArm] = deal(auxdata{:});

nA = M + nDof; % number of actuators
J = zeros(nDof*N,(nA)*N);
for i = 1:N
    for k = 1:nDof
        J(k+(i-1)*nDof,nA*(i-1)+1:nA*i) = MomentArm((i-1)*nA+1:i*nA,k)';
    end
end

J = sparse(J);
end

% ------------------------------------------------------------------
function J = jacobian_SO (x, auxdata)
[M, N, nDof, Fmax, Fpas, ID_data, MomentArm] = deal(auxdata{:});

nA = M + nDof; % number of actuators
J = zeros(nDof*N,(nA)*N);
for i = 1:N
    for k = 1:nDof
        J(k+(i-1)*nDof,nA*(i-1)+1:nA*i) = Fmax((i-1)*nA+1:i*nA)'.*MomentArm((i-1)*nA+1:i*nA,k)';
    end
end

J = sparse(J);
end


