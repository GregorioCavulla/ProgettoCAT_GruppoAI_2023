clc;
clear all;
close all;

% parametri progetto
rs=1.5;
rr=1.3;
K=200;
gamma=0.3;
beta=0.3;
alfa=0.2;
ms=0.7;
mr=0.15;
ns_e=100;
nr_e=100;

%% Ingresso, stato, uscita ed equilibrio

% poichè possiamo intervenire su cf(t) si ha u(t)=cf(t)  

% poichè nr(t) si ha y(t)=nr(t);

% fissiamo come stato x(t)=[ns(t),nr(t)]

% ricaviamo l'uscita d'equilibrio
y_e=nr_e;

% dunque ricaviamo lo stato di equilibrio
x_e=[ns_e,nr_e];

% ponendo f(t)=0 ricaviamo
u_e=0;

%% linearizzazione

% linearizzando il sistema attorno all'equilibrio si ottiene
A_e= [-1.05 -0.45; -0.35 -0.95];
B_e= [-90; 5];
C_e= [0 1];
D_e= 0;

% otteniamo dunque il seguente sistema linearizzato
sis_lin=ss(A_e,B_e,C_e,D_e);

%% Funzione di trasferimento
% ricaviamo la funzione di trasferimento del sistema linearizzato
GG=tf(sis_lin);

% stampiamo su terminale la funzione di trasferimento mettendo in evidenza
% guadagno zeri e poli
zpk(GG)

% eseguiamo il plot di G
figure(1);
bode(GG);

% eseguiamo il plot degli zeri e dei poli di G
figure(2);
pzmap(GG);

%% Mapping specifiche

% 1) errore a regime nullo con riferimento a gradino
e_reg=0;

% 2) margine fase >= 40 deg
Mf_spec=40;

% 3) sovraelungazione percentuale massima del 7%
S_100_spec=0.07;
% questa specifica eqivale a dire Mf>=64.61, che è più stringente rispetto
% a 2) dunque ridefiniamo Mf_spec
Mf_spec=64.61;

% 3) tempo di assestamento al 5% inferiore a ad un secondo
T_a5_spec=1;

% 4) disturbo in uscita attenuato di almeno 60db in [0,0.05]
A_d=60;
omega_d_min=0;
omega_d_MAX=0.05;

% 5) disturbo in ingresso attenuato di almeno 90db in [10^4,10^6]
A_n=60;
omega_n_min=1e4;
omega_n_MAX=1e6;

% eseguiamo il plot della funzione di GG con le patch derivanti dalle
% specifiche
figure(3);
hold on;

% Specifiche su d
omega_d_min = 0.0001; % lower bound per il plot
Bnd_d_x = [omega_d_min; omega_d_MAX; omega_d_MAX; omega_d_min];
Bnd_d_y = [A_d; A_d; -150; -150];
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche su n
Bnd_n_x = [omega_n_min; omega_n_MAX; omega_n_MAX; omega_n_min];
Bnd_n_y = [-A_n; -A_n; 100; 100];
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche tempo d'assestamento
omega_Ta_min =  0.0001; % lower bound per il plot
omega_Ta_MAX = 300 /(Mf_spec*T_a5_spec);
Bnd_Ta_x = [omega_Ta_min; omega_Ta_MAX; omega_Ta_MAX; omega_Ta_min];
Bnd_Ta_y = [0; 0; -150; -150];
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_mag = ["A_d", "A_n", "\omega_{c,min}"];
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(GG);
grid on; zoom on;

% Specifiche sovraelongazione (margine di fase)
omega_c_min = omega_Ta_MAX;
omega_c_MAX = omega_n_min;

phi_spec = Mf_spec - 180;
phi_low = -270; % lower bound per il plot

Bnd_Mf_x = [omega_c_min; omega_c_MAX; omega_c_MAX; omega_c_min];
Bnd_Mf_y = [phi_spec; phi_spec; phi_low; phi_low];
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
hold on;

% Legenda colori
Legend_arg = ["G(j\omega)","M_f"];
legend(Legend_arg);

%% Sintesi regolatore statico
s= tf("s");
R_s=1/s; % scegliamo di non fissare mu_s

%% Sintesi sistema esteso
GG_e=R_s*GG;

% eseguiamo il plot della funzione di GG_e con le patch derivanti dalle
% specifiche
figure(4);
hold on;

% Specifiche su d
omega_d_min = 0.0001; % lower bound per il plot
Bnd_d_x = [omega_d_min; omega_d_MAX; omega_d_MAX; omega_d_min];
Bnd_d_y = [A_d; A_d; -150; -150];
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche su n
Bnd_n_x = [omega_n_min; omega_n_MAX; omega_n_MAX; omega_n_min];
Bnd_n_y = [-A_n; -A_n; 100; 100];
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche tempo d'assestamento
omega_Ta_min =  0.0001; % lower bound per il plot
omega_Ta_MAX = 300 /(Mf_spec*T_a5_spec);
Bnd_Ta_x = [omega_Ta_min; omega_Ta_MAX; omega_Ta_MAX; omega_Ta_min];
Bnd_Ta_y = [0; 0; -150; -150];
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_mag = ["A_d", "A_n", "\omega_{c,min}"];
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(GG_e);
grid on; zoom on;

% Specifiche sovraelongazione (margine di fase)
omega_c_min = omega_Ta_MAX;
omega_c_MAX = omega_n_min;

phi_spec = Mf_spec - 180;
phi_low = -270; % lower bound per il plot

Bnd_Mf_x = [omega_c_min; omega_c_MAX; omega_c_MAX; omega_c_min];
Bnd_Mf_y = [phi_spec; phi_spec; phi_low; phi_low];
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
hold on;

% Legenda colori
Legend_arg = ["G_e(j\omega)","M_f"];
legend(Legend_arg);

%% Sintesi regolatore dinamico

Mf_star = Mf_spec+20; 
omega_c_star = 100;

% utilizziamo le formule di inversione
[mag_omega_c_star, arg_omega_c_star, ~] = bode(GG_e, omega_c_star);
mag_omega_c_star_db = 20*log10(mag_omega_c_star);

M_star = 10^(-mag_omega_c_star_db/20);
phi_star = Mf_star - 180 - arg_omega_c_star;

alpha_tau = (cos(phi_star*pi/180) - 1/M_star)/(omega_c_star*sin(phi_star*pi/180));
tau = (M_star - cos(phi_star*pi/180))/(omega_c_star*sin(phi_star*pi/180));
alpha = alpha_tau / tau;

% inseriamo un polo ad alte frequenze per rispettare la specifica sul
% disturbo di misura
tau_p = 1/500;
polo_hf = 1/(1 + tau_p*s);

% sintesi regolatore dinamico attraverso una rete anticipatrice
R_d = (1 + tau*s)/(1 + alpha_tau*s)*polo_hf;

% controllo di fisica realizzabilità
check_flag = min(tau, alpha_tau);

if check_flag < 0
    disp('Errore: polo/zero positivo');
    return;
end

%% Diagrammi di Bode sistema in anello chiuso con regolatore

% regolatore
RR = R_s*R_d;
% funzione di anello
LL = R_d*GG_e;

figure(5);
hold on;

% specifiche su ampiezza
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_n_x, Bnd_n_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);

% legenda colori
Legend_mag = ["A_d", "A_n", "\omega_{c,min}"];
legend(Legend_mag);

% plot Bode con margini di stabilità
margin(LL);
grid on; zoom on;

% specifiche su fase
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% legenda colori
Legend_arg = ["L(j\omega)","M_f"];
legend(Legend_arg);

%% Test sul sistema linearizzato

% funzioni di sensitività
FF=LL/(1+LL);
SS=1/(1+LL);

% calcolo il periodo di campionamento considerando che la frequenza deve
% essere almeno doppia rispetto alla frequenza massima
t_simulazione_n=0.02;
passo_n=0.5e-6;
tt_n=[0:passo_n:t_simulazione_n];

t_simulazione_d=3e3;
passo_d=5;
tt_d=[0:passo_d:t_simulazione_d];

t_simulazione=3;
passo=passo_n;
tt=[0:passo:t_simulazione];

% costruisco il segnale di riferimento
WW=-2;
ww=WW*ones(length(tt),1);

% costruisco i disturbi
DD=0.3; dd=0; wd=0.0125; dd_plot=0;
NN=0.2; nn=0; wn=1e4;   nn_plot=0;

for k=[1:1:4]
    dd=dd + DD*sin(wd*k*tt); 
    nn=nn + NN*sin(wn*k*tt); 
    dd_plot=dd_plot + DD*sin(wd*k*tt_d);
    nn_plot=nn_plot+NN*sin(wn*k*tt_n);  
end 

% calcolo l'uscita per il segnale di riferimento
y_w=lsim(FF,ww,tt);

% calcolo le uscite per il segnale di riferimento e per i disturbi
y_d=lsim(SS, dd, tt);
y_d_plot=lsim(SS, dd_plot, tt_d);

y_n=lsim(-FF, nn, tt);
y_n_plot=lsim(-FF,nn_plot,tt_n);

% eseguo il plot per la risposta al segnale di riferimento
figure(6);
grid on, zoom on, hold on;
plot(tt, ww, 'm');
plot(tt, y_w, 'b');
legend('ww','y_w');

% eseguo il plot per la risposta al disturbo di attuazione
figure(7);
hold on, grid on, zoom on
plot(tt_d,dd_plot,'m');
plot(tt_d,y_d_plot,'b');
legend('dd','y_d');

% eseguo il plot per la risposta al disturbo di misura
figure(8);
hold on, grid on, zoom on;
plot(tt_n,nn_plot,'m');
plot(tt_n,y_n_plot,'b');
legend('nn','y_n');

% eseguo il plot sul sistema linearizzato considerando sia il segnale di
% riferimento che i rumori
y_tot = y_w + y_d + y_n;
figure(9);
grid on, zoom on, hold on;
plot(tt,ww,'m');
plot(tt, y_tot,'b');

patch([0,t_simulazione,t_simulazione,0],[WW*(1+S_100_spec),WW*(1+S_100_spec),WW-1,WW-1],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);
patch([T_a5_spec,t_simulazione,t_simulazione,T_a5_spec],[WW*(1-0.05),WW*(1-0.05), WW+1, WW+1],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_a5_spec,t_simulazione,t_simulazione,T_a5_spec],[WW*(1+0.05),WW*(1+0.05),WW-1, WW-1],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);
legend('ww','y_tot');









