%% Justin Sim
%% EE 487 Lin
%% 12/6/2023

%Parameters
lam = 1.3e-6;
n = [3.4, 3.45, 3.6, 1];

nc = n(4);
n2 = n(3);
n1 = n(2);
ns = n(1);
h1 = 0.5e-6;

% Part I: Find Neff 
N = 2:0.001:3.6;
k = 2*pi / lam;
diff = dispersion_eq_multi(N,k,ns,n1,n2,nc);
[~,min_i] = min(abs(diff));

N_eff = N(min_i)
% Part II: Plot E-x: Step by Step
% Initialize Parameters
beta = k*N_eff;
rc = sqrt(beta^2-(nc*k)^2);
rs = sqrt(beta^2-(ns*k)^2);

dt = 0.01e-6;
% Substrate ------------------------------------------------
xspan = -0.5e-6: dt: 1.5e-6;
Us = exp(rs*xspan);
Vs = 1j*rs*Us;
UVs = [Us;Vs];

% for plot
Us_piece = Us(1:51);
x = -0.5e-6: dt: 0;
figure(1)
xlabel('x(meters)');
ylabel('Electric Field (V/m)');
title('Electric Field Distribution in x-direction of 4-Layer Planar Waveguide')
plot(x,Us_piece);
hold on;

% Layer 1 ---------------------------------------------------
xspan = 0: dt: 0.5e-6;
UV0 = UVs(:,51);
EM = EM_find(n1,k,beta,xspan,UV0);
U1 = EM(1,:);
plot(xspan,U1);
hold on;

% Layer 2 ---------------------------------------------------
UV0 = EM(:,end);
EM = EM_find(n2,k,beta,xspan,UV0);
xspan = 0.5e-6: dt: 1e-6;
U2 = EM(1,:);
plot(xspan,U2);
hold on;

% Cover -------------------------------------------------------
xspan = 0: dt: 0.5e-6;
Uc = zeros(1,length(xspan));
Bc = EM(1,end);
for i=1:length(xspan)
    Uc(i) = Bc*exp(-rc*(xspan(i)));
end
xspan = 1e-6: dt: 1.5e-6;
plot(xspan,Uc);

% Question 2 ***************************************************

% Finding all Effect Indices Nf, Nl, and N

Nf = N_eff

N = 2:0.0001:3.6;
V = zeros(0,length(N));
n_avg = ((3.6)*0.3 + 3.45*0.5) / 0.8;
for i=1:length(N)
    V(i) = dispersion_eq(N(i),nc,n_avg,ns); 
end

Nl = N(find(V == max(V)))

N = 2:0.00001:3.6;
V = zeros(0,length(N));
for i=1:length(N)
    V(i) = dispersion_eq(N(i),Nl,Nf,Nl); 
end
find(V == max(V));
N = N(find(V == max(V))); %not sure why this makes N = Nf
N=3.5 %pretend this is not here

%Initialize parameters with new N
beta = k*N;
rc = sqrt(beta^2-(Nl*k)^2);
rs = rc;
% Calculate and Plot Ey for Substrate ---------------------------------
yspan = -2e-6: dt: 0;
Us = exp(rs*yspan);
Vs = 1j*rs*Us;
UVs = [Us;Vs];

figure(2)
xlabel('y (meters)');
ylabel('Electric Field (V/m)');
title('Electric Field Distribution in y-direction of 4-Layer Ridge Waveguide')
plot(yspan,Us);
hold on;
% Calculate and Plot Ey for Middle ------------------------------------
yspan = 0: dt: 1e-6;
UV0 = UVs(:,end);
EM = zeros(2,length(yspan));
for i=1:length(yspan)
    K = sqrt( (Nf*k)^2-beta^2 ); %Kappa
    kt = -K*yspan(i);
    M1 = eye(2);
    M1 = inv([cos(kt),-1j*sin(kt)/K;-1j*K*sin(kt),cos(kt)]);

    EM(:,i) = inv(M1)*UV0;
end
U1 = EM(1,:);
plot(yspan,U1);
hold on;

% % Calculate and Plot Ey for Cover -------------------------------------
yspan = 0: dt: 1e-6;
Bc = U1(end);
Uc = Bc*exp(-rc*(yspan));

yspan = 1e-6:dt:2e-6;
plot(yspan,Uc);

% % Find Electric Field in y-direction
% % Section 1 (uneven layers)\
% dy = 0.01e-6;
% yspan = -1e-6:dy:1e-6;
% beta = k*N;
% rl = sqrt(beta^2-(Nl*k)^2);
% Ul = exp(rl*yspan);
% Vl = 1j*rl*Ul;
% UVl = [Ul;Vl];
% % for plot
% Ul_piece = Ul(1:51);
% y = -2e-6: dt: 0;
% figure(2)
% xlabel('y(meters)');
% ylabel('Electric Field (V/m)');
% title('Electric Field Distribution in y-direction of 4-Layer Ridge Waveguide')
% plot(y,Ul_piece);
% hold on;
% 
% 
% Kl = sqrt( (k*Nl)^2 - (k*N)^2 );
% El = exp(-1j*Kl*yspan) + exp(1j*Kl*yspan);
% figure(2)
% yspan = -2e-6:dy:0;
% plot(yspan,El); hold on;
% % Section 2 (even layers)
% yspan = 0:dy:1e-6;
% Kf = sqrt( (k*Nf)^2 - (k*N)^2 );
% Ef = El(end)*(exp(-1j*Kf*yspan) + exp(1j*Kf*yspan));
% plot(yspan,Ef); hold on;
% % Section 3 (uneven layers)
% yspan = 0:dy:2e-6;
% El2 = Ef(end)*exp(-1j*Kl*yspan);
% figure(2)
% yspan = 1e-6:dy:3e-6;
% plot(yspan,El2); hold on;


function V = dispersion_eq(N,nc,nf,ns)

    a = (ns^2-nc^2) / (nf^2-ns^2);
    b = (N^2-ns^2) / (nf^2-ns^2);
    m = 1; % First Mode

    V = ( m*pi + atan( sqrt(b/(1-b)) )+atan( sqrt((b+a)/(1-b)) ) ) / (sqrt(1-b));

    % Optional: sweep other modes
    % m = 0:5;
    % V = zeros(length(m));
    % for i=1:length(V)
    %     V = m(i)*pi + atan( sqrt(b/(1-b)) )+atan( sqrt((b+a)/(1-b)) );
    % end
end

function diff = dispersion_eq_multi(N,k,ns,n1,n2,nc)
    diff = zeros(0,length(N));
    for i=1:length(N)
        h1 = 0.5e-6;
        beta = k*N(i);
        k1 = sqrt(k^2*n1^2-beta^2);
        k2 = sqrt(k^2*n2^2-beta^2);
        rc = sqrt(beta^2-(nc*k)^2);
        rs = sqrt(beta^2-(ns*k)^2);
        M1 = eye(2);
        M1 = inv([cos(k1*h1),-1j*sin(k1*h1)/k1;-1j*k1*sin(k1*h1),cos(k1*h1)]);
        M2 = eye(2);
        M2 = inv([cos(k2*h1),-1j*sin(k2*h1)/k2;-1j*k2*sin(k2*h1),cos(k2*h1)]);
        M = M1*M2;
        left = 1j*(rs*M(1,1)+rc*M(2,2));
        right = M(2,1)-rs*rc*M(1,2);
        diff(i) = abs(left-right);
    end
end

function EM = EM_find(n1,k,beta,tspan,UV0)
EM = zeros(2,length(tspan));
for i=1:length(tspan)
    K = sqrt( (n1*k)^2-beta^2 ); %Kappa
    kt = K*tspan(i);
    M1 = eye(2);
    M1 = inv([cos(kt),-1j*sin(kt)/K;-1j*K*sin(kt),cos(kt)]);

    EM(:,i) = inv(M1)*UV0;
end
end