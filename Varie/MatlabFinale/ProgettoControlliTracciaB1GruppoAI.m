clc; 
clear all; 
close all;

% Dichiarazione delle variabili globali
global Rs Rr K Ms Mr beta alfa gamma;

% Definizione delle costanti
Rs = 1.5;       % tasso di riproduzione delle cellule suscettibili
Rr = 1.3;       % tasso di riproduzione delle cellule resistenti
Ms = 0.7;       % tasso di mortalità delle cellule suscettibili
Mr = 0.15;      % tasso di mortalità delle cellule resistenti
K = 200;        % massimo numero di cellule che possono esistere all'interno dell'ambiente
gamma = 0.3;    % termine di mutazione r -> s
beta = 0.3;     % termine di mutazione s -> r
alfa = 0.2;     % termine di mutazione s -> r causata dal trattamento terapeutico

% Definizione delle condizioni iniziali
x_1e = 100;     % popolazione iniziale di cellule suscettibili
x_2e = 100;     % popolazione iniziale di cellule resistenti
x0 = [x_1e; x_2e];  % vettore contenente le popolazioni iniziali delle due tipologie di cellule

%%Linearizzazione

% Troviamo il valore di equilibrio per l'ingresso u_e ponendo a 0 la f(x) e troviamo:
u_e = 0;
y_e = x_2e;

% Matrici del sistema linearizzato
A_eq = [-1.05 -0.45; -0.35 -0.95];
B_eq = [-90; 5];
C_eq = [0 1];
D_eq = 0;

% Costruiamo il modello nello spazio degli stati
modello = ss(A_eq, B_eq, C_eq, D_eq);

%% punto 2: funzione di trasferimento

% Costruzione della funzione di trasferimento G a partire dal modello di stato
G = tf(modello);

% Zeri, poli e guadagno della funzione di trasferimento
zpk(G)


figure(1);
bode(G); % Plot del diagramma di Bode


figure(2);
pzmap(G); % Plot zeri e poli

%% Punto 3: definizione specifiche

% ampiezze 
global DD NN;

WW = -2;    %gradino  
DD = 0.3;   %sinusoide
NN = 0.2;   %sinusoide

% errore a regime
e_star = 0;     % sarà necessario un polo nell'origine

% attenuazione disturbo sull'uscita
A_d = 60;
omega_d_min = 0.0001; % lower bound per il plot
omega_d_max = 0.05;

% attenuazione disturbo di misura
A_n = 90;
omega_n_min = 1e4;
omega_n_max = 1e6;

% Sovraelongazione Massima
S_100_spec = 0.07;

% Margine di fase
logsq = (log(S_100_spec))^2;
xi = sqrt(logsq/(pi^2+logsq));
Mf_spec = xi*100;

% Tempo di assestamento al 5%
T_a5_spec = 1;

% limite inferiore per wc dato da specifica sul Ta
omega_Ta_min = 1e-4; % lower bound per il plot
omega_Ta_max = 300/(Mf_spec*T_a5_spec); 

%% Regolatore statico e sistema esteso Ge

s = tf('s');
mu_s = 1; % lo metto NON per errore a regime ma per specifica sul Ta

R_s = mu_s/s;   % polo nell'origine per errore a regime
Ge = G*R_s;

%% patch e diagramma di bode di Ge

figure(3);
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

% label colori
label_mag = ["A_d", "A_n", "\omega_{c,min}"];
legend(label_mag);

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

% label colori
label_arg = ["G_e(j\omega)","M_f"];
legend(label_arg);

%% Design del regolatore dinamico

Mf_star = Mf_spec + 12;
omega_c_star = 100;

% formule di inversione (anticipatore)
[mag_omega_c_star, arg_omega_c_star, ~] = bode(Ge, omega_c_star);
arg_omega_c_star;
mag_omega_c_star_db = 20*log10(mag_omega_c_star);

M_star = 10^(-mag_omega_c_star_db/20);
phi_star = Mf_star - 180 - arg_omega_c_star;
phi_star_rad = deg2rad(phi_star);

alpha_tau =  (cos(phi_star_rad) -1/M_star)/(omega_c_star*sin(phi_star_rad));
tau = (M_star - cos(phi_star_rad))/(omega_c_star*sin(phi_star_rad));
alpha = alpha_tau / tau;

% controllo realizzabilità della rete anticipatrice
if min(tau, alpha_tau) < 0
    disp('Errore: polo/zero positivo');
    return;
end

%% Diagrammi di Bode con specifiche includendo regolatore dinamico
tau_p = 1/500;
polo_hf = 1/(1 + tau_p*s);  % polo in alta frequenza per rispettare la specifica su n(t)

R_d = (1 + tau*s)/(1 + alpha_tau*s)*polo_hf; % rete anticipatrice
R = R_s*R_d; % regolatore
LL = R_d*Ge; % funzione di anello

figure(4);
hold on;

% Specifiche su ampiezza
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_n_x, Bnd_n_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);

% label colori
label_mag = ["A_d", "A_n", "\omega_{c,min}"];
legend(label_mag);

% Plot Bode con margini di stabilità
margin(LL);
grid on; zoom on;

% Specifiche su fase
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% label colori
label_arg = ["L(j\omega)"; "M_f"];
legend(label_arg);


%% Punto 4: testare il modello linearizzato

% funzioni di sensitività
FF= LL/(1+LL);
SS=1/(1+LL);

%% risposta al riferimento y_w

T_simulazione = 3;
passo = 0.7e-6; % sample rate doppio rispetto alla frequenza massima 

global tt ww;
tt = (0:passo:T_simulazione);

% costruisco il riferimento
ww = WW * ones(length(tt), 1);
y_w = lsim(FF, ww, tt);

figure(5);
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
% tt_n = (0:0.7e-6:1e-2);  

for k=[1:1:4]
    dd = dd + DD*sin(wd*k*tt); 
    nn = nn + NN*sin(wn*k*tt); 
end

% uscite dei rumori
y_d = lsim(SS, dd, tt);
y_n = lsim(-FF, nn, tt);

% plot per d
figure(6);
hold on, grid on, zoom on
plot(tt,dd,'m')
plot(tt,y_d,'b')
legend('dd','y_d')

% plot per n
figure(7);
hold on, grid on, zoom on
plot(tt,nn,'m')
plot(tt,y_n,'b')
legend('nn','y_n')

% plot risposta al gradino + rumori 
y_tot = y_w + y_d + y_n;
figure(8);
plot(tt, y_tot,'b')

% vincolo sovraelongazione
patch([0,T_simulazione,T_simulazione,0],[WW*(1+S_100_spec),WW*(1+S_100_spec),WW-1,WW-1],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);

% vincolo tempo di assestamento all'5%
patch([T_a5_spec,T_simulazione,T_simulazione,T_a5_spec],[WW*(1-0.05),WW*(1-0.05), WW+1, WW+1],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_a5_spec,T_simulazione,T_simulazione,T_a5_spec],[WW*(1+0.05),WW*(1+0.05),WW-1, WW-1],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

label_step = ["y_{tot}"; "Vincolo sovraelongazione"; "Vincolo tempo di assestamento"];
legend(label_step);
