clc; clear all; close all;

global Rs Rr K Ms Mr beta alfa gamma y_e WW u_e nn dd;

Rs=1.5; Rr=1.3;     %tassi di riproduzione
Ms=0.7; Mr=0.15;    %tassi di mortalità
K=200;              %max cellule dentro ambiente
gamma=0.3;          %termine di mutazione r->s
beta=0.3;           %termine di mutazione s->r
alfa=0.2;           %termine di mutazione s->r causata dal trattamento

x_1e=100; x_2e=100;   %coppia di equilibrio
x0 = [x_1e; x_2e];    %stato iniziale
interv = [0 100];    

Cf = @(t) 1;        %ingresso di controllo u(t)

%% punto 1: forma di stato
dyn = @(t, N) [-Rs*log((N(1)+N(2))/K)*N(1)-(Ms*Cf(t)*N(1))-(beta*N(1))+(gamma*N(2))-(alfa*Cf(t)*N(1));
               -Rr*log((N(1)+N(2))/K)*N(2)-(Mr*Cf(t)*N(2))+(beta*N(1))-(gamma*N(2))+(alfa*Cf(t)*N(1))];

[time, traj] = ode45(dyn, interv, x0);

yy=traj(:,2);   %y=x(2)

%plot dell'evoluzione nel tempo
if 0
    figure(1);
    plot(time, traj);
    %plot(time, yy);  %possiamo misurare solo quelle resistenti;
    grid on; zoom on;
    legend("cellule suscettibili", "cellule resistenti");
    %legend("cellule resistenti");
end

%% Linearizzazione

%cerchiamo Ue ingresso di equilibrio ponendo a 0 la f(x) e troviamo
u_e=0;
y_e=x_2e;

%matrici del sistema linearizzato
A_eq = [-1.05 -0.45; -0.35 -0.95];
B_eq = [-90; 5];
C_eq = [0 1];
D_eq = 0;

% state-space model
modello = ss(A_eq, B_eq, C_eq, D_eq);

% plot del sistema linearizzato
if 0
    cf_vettore = ones(length(interv), 1);
    [YY, TT, XX] = lsim(modello, cf_vettore, interv, x0);
    
    figure(2);
    plot(TT, XX);  
    grid on; zoom on;
    %legend("cellule resistenti");
    legend("cellule suscettibili", "cellule resistenti");
end

%% punto 2: funzione di trasfermento

G = tf(modello);
zpk(G)
    
%diagramma di bode e pz-map di G
if 0
    figure(3);
    bode(G);
    figure(4);
    pzmap(G);
end

%% Punto 3: definizione specifiche

% ampiezze 
global DD NN;

WW = -2;    %gradino  
DD = 0.3;   %sinusoide
NN = 0.2;   %sinusoide

% errore a regime
e_star = 0;     % abbiamo bisogno di un polo nell'origine

% attenuazione disturbo sull'uscita
A_d = 60;
omega_d_min = 0.0001; % lower bound per il plot
omega_d_max = 0.05;

% attenuazione disturbo di misura
A_n = 90;
omega_n_min = 1e4;
omega_n_max = 1e6;

% Sovraelongazione massima e tempo d'assestamento all'1%
S_100_spec = 0.07;

% Margine di fase
logsq = (log(S_100_spec))^2;
xi = sqrt(logsq/(pi^2+logsq));
Mf_spec = xi*100

% Tempo di assestamento
T_a5_spec = 1;

% limite inferiore per wc dato da specifica sul Ta
omega_Ta_min = 1e-4; % lower bound per il plot
omega_Ta_max = 300/(Mf_spec*T_a5_spec); 

%% Regolatore statico e sistema esteso Ge

s = tf('s');
mu_s = 1; % lo metto NON per errore a regime ma per specifica sul Ta

R_s = mu_s/s;   % la G NON ha un polo nell'origine quindi lo dobbiamo mettere qua!
Ge = G*R_s;

%% patch e diagramma di bode di Ge

figure(5);
hold on;

% Specifiche su d
Bnd_d_x = [omega_d_min; omega_d_max; omega_d_max; omega_d_min];
Bnd_d_y = [A_d; A_d; -150; -150];
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche su n
Bnd_n_x = [omega_n_min; omega_n_max; omega_n_max; omega_n_min];
Bnd_n_y = [-A_n; -A_n; 100; 100];
patch(Bnd_n_x, Bnd_n_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche tempo d'assestamento (minima pulsazione di taglio)
Bnd_Ta_x = [omega_Ta_min; omega_Ta_max; omega_Ta_max; omega_Ta_min];
Bnd_Ta_y = [0; 0; -150; -150];
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_mag = ["A_d", "A_n", "\omega_{c,min}", "G_e(j\omega)"];
legend(Legend_mag);

bode(Ge);

grid on; zoom on;

% Specifiche sovraelongazione (margine di fase)
omega_c_min = omega_Ta_max;
omega_c_max = omega_n_min;

phi_spec = Mf_spec - 180;
phi_low = -270; % lower bound per il plot

Bnd_Mf_x = [omega_c_min; omega_c_max; omega_c_max; omega_c_min];
Bnd_Mf_y = [phi_spec; phi_spec; phi_low; phi_low];
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_arg = ["G_e(j\omega)"; "M_f"];
legend(Legend_arg);

%% Design del regolatore dinamico

% qua c'è del trial & error
Mf_star = Mf_spec+10;
omega_c_star = 100;

% scrivere formule di inversione (anticipatore)
[mag_omega_c_star, arg_omega_c_star, ~] = bode(Ge, omega_c_star);
arg_omega_c_star
mag_omega_c_star_db = 20*log10(mag_omega_c_star)

M_star = 10^(-mag_omega_c_star_db/20)
phi_star = Mf_star - 180 - arg_omega_c_star;
phi_star_rad = deg2rad(phi_star)

alpha_tau =  (cos(phi_star_rad) -1/M_star)/(omega_c_star*sin(phi_star_rad))
tau = (M_star - cos(phi_star_rad))/(omega_c_star*sin(phi_star_rad))
alpha = alpha_tau / tau

%controllo realizzabilità della rete anticipatrice
check_flag = min(tau, alpha_tau);
if check_flag < 0
    disp('Errore: polo/zero positivo');
    return;
end

%% Diagrammi di Bode con specifiche includendo regolatore dinamico
tau_p = 1/600;
polo_hf = 1/(1 + tau_p*s);  % polo in alta frequenza per rispettare la specifica su n(t)

R_d = (1 + tau*s)/(1 + alpha_tau*s)*polo_hf; % rete anticipatrice
LL = R_d*Ge; % funzione di anello

figure(6);
hold on;

% Specifiche su ampiezza
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_n_x, Bnd_n_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_mag = ["A_d", "A_n", "\omega_{c,min}", "L(j\omega)"];
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(LL);
grid on; zoom on;

% Specifiche su fase
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_arg = ["L(j\omega)"; "M_f"];
legend(Legend_arg);

%% Punto 4: testare il modello linearizzato

% funzioni di sensitività
FF= LL/(1+LL);
SS=1/(1+LL);

%% risposta al riferimento y_w

T_simulazione = 3;
passo = 0.5e-6; % devo campionare al doppio della Fmax del segnale a frequenza maggiore

global tt ww;
tt = (0:passo:T_simulazione);

% costruisco il riferimento
ww = WW * ones(length(tt), 1);
y_w = lsim(FF, ww, tt);

figure(7);
grid on, zoom on, hold on;

plot(tt, ww, 'm');
plot(tt, y_w, 'b');
legend('ww','y_w')

%% risposte ai rumori

% costruisco i rumori 
global wd wn;

dd=0;   wd=0.0125;
nn=0;   wn=1e4;

% tt_d = (0:1e-2:1e3);
% tt_n = (0:0.5e-6:1e-2);  

for k=[1:1:4]
    dd = dd + DD*sin(wd*k*tt); 
    nn = nn + NN*sin(wn*k*tt); 
end

% uscite dei rumori
y_d = lsim(SS, dd, tt);
y_n = lsim(-FF, nn, tt);

%plot per d
figure(8);
hold on, grid on, zoom on
plot(tt,dd,'m')
plot(tt,y_d,'b')
legend('dd','y_d')

%plot per n
figure(9);
hold on, grid on, zoom on
plot(tt,nn,'m')
plot(tt,y_n,'b')
legend('nn','y_n')

%plot risposta al gradino + rumori 
y_tot = y_w + y_d + y_n;
figure(10);
plot(tt, y_tot,'b')

% vincolo sovraelongazione
patch([0,T_simulazione,T_simulazione,0],[WW*(1+S_100_spec),WW*(1+S_100_spec),WW-1,WW-1],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);

% vincolo tempo di assestamento all'5%
LV = WW; % in questo caso non c'è errore a regime quindi possiamo fare così
patch([T_a5_spec,T_simulazione,T_simulazione,T_a5_spec],[LV*(1-0.05),LV*(1-0.05), LV+1, LV+1],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_a5_spec,T_simulazione,T_simulazione,T_a5_spec],[LV*(1+0.05),LV*(1+0.05),LV-1, LV-1],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

Legend_step = ["y_{tot}"; "Vincolo sovraelongazione"; "Vincolo tempo di assestamento"];
legend(Legend_step);

%% Punto 5: test sul sistema NON lineare

if 1
    global A_R B_R C_R D_R;
    
    R = R_s*R_d;
    [A_R, B_R, C_R, D_R] = tf2ss(cell2mat(R.num), cell2mat(R.den));
    
    x0 = [100; 100; 0; 0; 0];  % ridefinisco lo stato iniziale che adesso deve essere 4-dimensionale
    [time_full, traj_full] = ode45(@sistemaNonLineareChiuso, tt, x0);
    
    % uscita più disturbo in uscita
    y = traj_full(:, 2);
    
    figure(11)
    hold on; zoom on; grid on;
    plot(time_full, y, 'b');

    legend('y_w nonLineare')
end 

function out = sistemaNonLineareChiuso(t, x)
    % devo ridefinire anche qua dentro le var globali se no warning
    global Rs Rr K Ms Mr beta alfa gamma y_e WW u_e A_R B_R C_R D_R tt DD NN wd wn;
    
    dd=0; nn=0;

    % mi ricostruisco i rumori con periodo di ode45 t
    for k=[1:1:4]
        dd = dd + DD*sin(wd*k*t); 
        nn = nn + NN*sin(wn*k*t); 
    end 

    % output del sistema y(t) con disurbo in uscita
    y = x(2) + dd;
    % input del regolatore e(t) con disturbo di misura
    e = y_e + WW - y - nn;
    %input del sistema, ovvero output del regolatore u(t)
    u = (C_R*x(3:5, 1) + D_R*e) + u_e;


    % out contiene sia gli stati del sistema che gli stati del regolatore
    % in modo da far calcolare ad ode45 entrambi
    out = [-Rs*log((x(1)+x(2))/K)*x(1)-(Ms*u*x(1))-(beta*x(1))+(gamma*x(2))-(alfa*u*x(1));
           -Rr*log((x(1)+x(2))/K)*x(2)-(Mr*u*x(2))+(beta*x(1))-(gamma*x(2))+(alfa*u*x(1));
           A_R(1, :)*x(3:5, 1)+B_R(1, :)*e;
           A_R(2, :)*x(3:5, 1)+B_R(2, :)*e;
           A_R(3, :)*x(3:5, 1)+B_R(3, :)*e];
end

