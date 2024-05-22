%--------------- Geometric constants --------------------------------------

l_b = 0.7;          % length of the beam
b = 0.014;          % width of the beam
h_b = 0.014;        % thickness of the beam

l_p = 0.067;        % length of a portion with piezo
h_p = 0.002;        % thickness of a piezo

I_b=b*(h_b^3)/12;   % second moment of area of the structure 
I_p=b*((h_b+2*h_p)^3-h_b^3)/12/2;  % second moment of area for one piezo

%--------------- Material constants ---------------------------------------

rho_b = 7850;             % density of the main structure  
Y_b = 210*10^9;           % Young modulus of the main structure 

rho_p = 7800;             % density of the piezo 
Y_p_E = 66e9;             % elastic stiffness of the piezo at constant electric field

% ==============================
n = 10;             % number of unit cells
a = l_b/n;          % length of a unit cell 

YI_E =a/(l_p/(Y_b*I_b+2*Y_p_E*I_p)+(a-l_p)/(Y_b*I_b));   % homogenized Young modulus*I 
rhoS = (rho_b*b*h_b*a+2*rho_p*b*h_p*l_p)/a;              % homogenized linear density rho*S

% ==============================



Le=a; % lenght of the element

Km=...
YI_E/(Le^3)*[12,   6*Le,   -12,   6*Le;...
            6*Le, 4*Le^2, -6*Le, 2*Le^2;...
            -12,  -6*Le, 12,  -6*Le;...
            6*Le, 2*Le^2, -6*Le, 4*Le^2];      % mechanical stiffness matrix      


Mm=...
rhoS*Le/420*[ 156,   22*Le, 54,  -13*Le;...
            22*Le,  4*Le^2,   13*Le, -3*Le^2;...
            54,   13*Le,    156,  -22*Le;...
            -13*Le, -3*Le^2,  -22*Le,  4*Le^2]; % mechanical mass matrix
                

% ==============================



%--------------- Full FEM matrices ------------------------------------------

n_dof_uc = 4;                     % number of degrees of freedom for a unit cell
n_dof_beam = (n_dof_uc-2)*n+2;    % number of degrees of freedom for the beam

K_full= zeros(n_dof_beam);        % start with empty matrices 
M_full= zeros(n_dof_beam);        

for id=1:n                        % fill the matrices depending on position of the element
    pos = (n_dof_uc-2)*(id-1)+1:(n_dof_uc-2)*id+2;             
    K_full(pos,pos)=K_full(pos,pos)+Km;
    M_full(pos,pos)=M_full(pos,pos)+Mm;
end

%--------------- Boundary conditions ------------------------------------------
   
% for a cantilever beam  
dofs_fixed = [1,2]';             % first displacement and rotation blocked

dofs_free=setdiff(1:n_dof_beam,dofs_fixed);   % keep only the free dof
K_tot = K_full(dofs_free,dofs_free);   
M_tot = M_full(dofs_free,dofs_free);


%--------------- Frequencies and mode shapes ------------------------------------------

[V,D] = eig(K_tot,M_tot);                     % compute the eigenproperties
[sortedD, I] = sort(diag(D)); 
freq_123 = sqrt(sortedD(1:3))/2/pi          % first three frequencies

disp_dof = 1:2:n_dof_beam-2;                  % degrees of freedom for displacement
figure
plot([0:10],[0;V(disp_dof,find(I == 1))])                % plot mode shapes
hold on
plot([0:10],[0;V(disp_dof,find(I == 2))])
plot([0:10],[0;V(disp_dof,find(I == 3))])
xlabel('Position');
ylabel('Mode shape');
box on



% ==============================


%--------------- External excitation ---------------------------------

F_ext = zeros(size(K_tot,2),1);     % external force vector
F_ext(size(K_tot,2)-1)=1;           % unitary force applied at the free end of the beam


%--------------- Parameters for the frequency response functions ------------------

f_min = 1;                      % minimum frequency
f_max = 500;                    % maximum frequency
nb_pt = 1000;                   % number of frequency points

om_min = f_min*2*pi;            % minimum angular frequency
om_max = f_max*2*pi;            % maximum angular frequency
x = linspace(om_min,om_max,nb_pt);   % angular frequency vector


%--------------- Frequency loops ---------------------------------

disp_free_end = zeros(size(x,2),1);     % displacement at the free end of the beam 
 
for j=1:size(x,2)
    om = x(j);                                   % angular frequency
    q_dof = (K_tot-om^2*M_tot)\F_ext;            % solve for forced excitation
    disp_free_end(j) =  q_dof(size(K_tot,2)-1);  % select degree of freedom
end
figure
plot(x/(2*pi),20*log10(abs(disp_free_end)),'-b','Linewidth',2)  % plot frequency responce function
hold on
xlabel('Frequency (Hz)');
ylabel('Compliance (dB, m/N)');
box on



%% Electrical beam 


C = 50.8e-9;    % capacitance
L = 0.247;      % inductance
a_tild = 2;     % transformer ratio 
R_L= 0;         % resistance in series with inductance
R_C= 0;         % resistance in series with capacitance
R_T = 838;      % resistance in transformers
n = 10;         % number of unit cells in the network


% ==============================



Cf = C*1e-3;   % numerical parameter for stiffness matrix

Ke=...       % electrical stiffness matrix
    4/a_tild^2/Cf*...
   [        1,                       a_tild/2,        -1,                       a_tild/2;...
     a_tild/2, (a_tild^2/4)*((C+2*Cf)/(C+Cf)), -a_tild/2,        (a_tild^2/4)*C/(C + Cf);...
           -1,                      -a_tild/2,         1,                      -a_tild/2;...
     a_tild/2,        (a_tild^2/4)*C/(C + Cf), -a_tild/2, (a_tild^2/4)*((C+2*Cf)/(C+Cf))];  
    

Me=...       % electrical mass matrix
    L/2*...   
    [        1 ,        a_tild/2,           0,              0;...
      a_tild/2 ,     a_tild^2/4 ,           0,              0;...  
              0,               0,          1 ,      -a_tild/2;...
              0,               0,  -a_tild/2 ,     a_tild^2/4];   
                                    
   
Ce=...       % electrical damping matrix;
    R_L/2*...   
    [        1 ,        a_tild/2,           0,              0;...
      a_tild/2 ,     a_tild^2/4 ,           0,              0;...  
              0,               0,          1 ,      -a_tild/2;...
              0,               0,  -a_tild/2 ,     a_tild^2/4]...
              +...      
    R_T/2*...   
    [        0 ,               0,           0,              0;...
             0 ,               1,           0,              0;...  
             0,               0,           0,              0;...
             0,               0,           0,              1]...
         +...      
    R_C*...   
    [        0 ,               0,           0,              0;...
             0 ,               1,           0,              -1;...  
             0,               0,           0,              0;...
             0,               -1,           0,              1];        
         



% ==============================



%%% (same code as for the mechanical beam with additional damping matrix) %%%

%--------------- Full FEM matrices ------------------------------------------

n_dof_uc = 4;                     % number of degrees of freedom for a unit cell
n_dof_beam = (n_dof_uc-2)*n+2;    % number of degrees of freedom for the beam

K_full= zeros(n_dof_beam);        % start with empty matrices 
C_full= zeros(n_dof_beam);
M_full= zeros(n_dof_beam);        

for id=1:n                        % fill the matrices depending on position of the element
    pos = (n_dof_uc-2)*(id-1)+1:(n_dof_uc-2)*id+2;             
    K_full(pos,pos)=K_full(pos,pos)+Ke;
    C_full(pos,pos)=C_full(pos,pos)+Ce;
    M_full(pos,pos)=M_full(pos,pos)+Me;
end

%--------------- Boundary conditions ------------------------------------------
   
% for a cantilever beam  
dofs_fixed = [1,2]';             % first displacement and rotation blocked

dofs_free=setdiff(1:n_dof_beam,dofs_fixed);   % keep only the free dof
K_tot = K_full(dofs_free,dofs_free); 
C_tot = C_full(dofs_free,dofs_free);  
M_tot = M_full(dofs_free,dofs_free);


%--------------- Frequencies and mode shapes ------------------------------------------

[V,D] = eig(K_tot,M_tot);                     % compute the eigenproperties
[sortedD, I] = sort(diag(D));                      
freq_123 = sqrt(sortedD(1:3))/2/pi            % first three frequencies

disp_dof = 1:2:n_dof_beam-2;                  % degrees of freedom for displacement
figure
plot([0:10],[0;V(disp_dof,find(I == 1))])                % plot mode shapes
hold on
plot([0:10],[0;V(disp_dof,find(I == 2))])
plot([0:10],[0;V(disp_dof,find(I == 3))])
xlabel('Position');
ylabel('Mode shape');
box on



% ==============================



%%% (same code as for the mechanical beam with additional damping matrix) %%%

%--------------- External excitation ---------------------------------

F_ext = zeros(size(K_tot,2),1);     % external force vector
F_ext(size(K_tot,2)-1)=1;           % unitary force applied at the free end of the beam


%--------------- Parameters for the frequency response functions ------------------

f_min = 1;                      % minimum frequency
f_max = 500;                    % maximum frequency
nb_pt = 1000;                   % number of frequency points

om_min = f_min*2*pi;            % minimum angular frequency
om_max = f_max*2*pi;            % maximum angular frequency
x = linspace(om_min,om_max,nb_pt);   % angular frequency vector


%--------------- Frequency loops ---------------------------------

disp_free_end = zeros(size(x,2),1);     % displacement at the free end of the beam 
 
for j=1:size(x,2)
    om = x(j);                                    % angular frequency
    q_dof = (K_tot+1i*om*C_tot-om^2*M_tot)\F_ext;            % solve for forced excitation
    disp_free_end(j) =  q_dof(size(K_tot,2)-1);  % select degree of freedom
end
figure
plot(x/(2*pi),20*log10(abs(disp_free_end)),'-b','Linewidth',2)  % plot frequency responce function
hold on
xlabel('Frequency (Hz)');
ylabel('Compliance (dB, m/N)');
box on

%% == COUPLING

e_th = 5.3e-3 ;        % global coupling coefficient
Kc=e_th*[0 1 0 -1]';   % coupling matrix

ndof = 4;              % number of degrees of freedom in the mechanical unit cell

K_coupled = [   Km+1/C*(Kc*Kc'),     1/C*Kc*[0, 1, 0, -1];...   % electromechanical stiffness matrix 
           1/C*[0;1;0;-1]*Kc',                        Ke];  
M_coupled = [Mm,zeros(ndof,4);zeros(4,ndof),Me];                     % electromechanical mass matrix 
C_coupled = [zeros(ndof,4),zeros(ndof,4);zeros(4,ndof),Ce];          % electromechanical damping matrix 
    

% ==============================

perm = zeros(ndof+4,ndof+4);     % build the permutation matrix        
perm(3:ndof,5:(ndof+2))=eye(ndof-2);   
perm(1,1)=1; perm(2,2)=1; perm(ndof+1,3)=1; perm(ndof+2,4)=1; perm(ndof+3,ndof+3)=1; perm(ndof+4,ndof+4)=1;  % permutation matrix to obtain the rearanged dynamic stiffness matrix
        
K_coupled=perm\K_coupled*perm;   % permutations to gather mechanical and electrical dofs
M_coupled=perm\M_coupled*perm;
C_coupled=perm\C_coupled*perm;


% ==============================


%%% (same code as for the mechanical beam with doubled degrees of freedom) %%%


%--------------- Full FEM matrices ------------------------------------------

n = 10;                           % number of unit cells 
n_dof_uc = 8;                     % number of degrees of freedom for a unit cell
n_dof_beam = (n_dof_uc-4)*n+4;    % number of degrees of freedom for the beam

K_full= zeros(n_dof_beam);        % start with empty matrices 
C_full= zeros(n_dof_beam);
M_full= zeros(n_dof_beam);        

for id=1:n                        % fill the matrices depending on position of the element
    pos = (n_dof_uc-4)*(id-1)+1:(n_dof_uc-4)*id+4;             
    K_full(pos,pos)=K_full(pos,pos)+K_coupled;
    C_full(pos,pos)=C_full(pos,pos)+C_coupled;
    M_full(pos,pos)=M_full(pos,pos)+M_coupled;
end

%--------------- Boundary conditions ------------------------------------------
   
% for a cantilever beam  
dofs_fixed = [1:4]';             % first displacement and rotation blocked

dofs_free=setdiff(1:n_dof_beam,dofs_fixed);   % keep only the free dof
K_tot = K_full(dofs_free,dofs_free); 
C_tot = C_full(dofs_free,dofs_free);  
M_tot = M_full(dofs_free,dofs_free);  


%--------------- External excitation ---------------------------------

F_ext = zeros(size(K_tot,2),1);     % external force vector
F_ext(size(K_tot,2)-3)=1;           % unitary force applied at the free end of the beam


%--------------- Parameters for the frequency response functions ------------------

f_min = 1;                      % minimum frequency
f_max = 500;                    % maximum frequency
nb_pt = 1000;                   % number of frequency points

om_min = f_min*2*pi;            % minimum angular frequency
om_max = f_max*2*pi;            % maximum angular frequency
x = linspace(om_min,om_max,nb_pt);   % angular frequency vector


%--------------- Frequency loops ---------------------------------

disp_free_end = zeros(size(x,2),1);     % displacement at the free end of the beam 
 
for j=1:size(x,2)
    om = x(j);                                    % angular frequency
    q_dof = (K_tot+1i*om*C_tot-om^2*M_tot)\F_ext;            % solve for forced excitation
    disp_free_end(j) =  q_dof(size(K_tot,2)-3);  % select degree of freedom
end
figure
plot(x/(2*pi),20*log10(abs(disp_free_end)),'-b','Linewidth',2)  % plot frequency responce function
hold on
xlabel('Frequency (Hz)');
ylabel('Compliance (dB, m/N)');
box on

% ==============================
