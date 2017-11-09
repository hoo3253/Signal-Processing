clc
clear all
close all

Folder='C:\Users\Hoo\Desktop\Abstract Animal Behavior\Matlab\Selected files - segundo análisis\';
File='C2_JCG3';

addpath(genpath('C:\Users\Hoo\Desktop\Abstract Animal Behavior\Matlab'));
[FileName,PathName,FilterIndex] = uigetfile('*.wav','Select a file',strcat(Folder,File));
addpath(genpath(PathName));
s_name=strcat(PathName,FileName);
n_fft=2^14;
window_length=0.05; %% en segundos
hop=50;
time_threshold=30;
pow_threshold=30;
f_i=350;
f_s=900;
rel_threshold=0.02;
pow_min=30;
%Perturbation measurements
ppq_num=5;% only odd numbers win_size

[s,fs]=audioread(s_name); %Audio reading
 s=s./abs(max(s));
% Parametros del modelo
 win_size=round(fs*window_length);% La menor frequencia de aleteo de aedes es 300 Hz, por lo tanto 0.014 s asegura al menos 4 ciclos de aleteo para eralizar la fft
 hop=round(win_size*hop/100);%Puntos de sobrelapamiento en porcentaje  
 w=hamming(win_size); %Ventana elegida para la transformada de fourier
 j=1;
 
 for t=1:hop:(length(s)-(n_fft))
    % Frecuencia
%      App(j)=max(s(t:(t+win_size-1)))-min(s(t:(t+win_size-1)));% Amplitud pico pico
     S=fft(w.*s(t:(t+win_size-1)),n_fft); % Transformada rápida de Fourier
     mod_S(:,j)=abs(S(1:length(S)/2,1));% Magnitud del espectro de potencia en la parte positiva
%      phase_S(:,j)=angle(S(1:length(S)/2,1));
%      phase_S_Unwraped(:,j)=unwrap(angle(S(1:length(S)/2,1)));
     dBFS(:,j)=20*log10(mod_S(:,j)/win_size); % dBFS Full scale, el valor maximo de la FFT sería el tamaño de la ventana usado. Este resultado muestra la señal en referencia a los puntos de la ventana usada. Por lo tanto el maximo 0dB sería tamaño de la ventana y de ahi para abajo - tantos dB
     j=j+1;
     
     t;
 end 
 
 f_inf=round(n_fft*f_i/fs);
 f_sup=round(n_fft*f_s/fs);
 
 %%% Vectores para plot frec
f=fs*linspace(0,1,n_fft);
f=f(1:length(f)/2);
 
for i=1:size(dBFS,2)        
    % Localización Frecuencias fundamentales 
    max_dBFS=max(dBFS(:,i));
    [pks_fund{i},locs]=findpeaks(dBFS(f_inf:f_sup,i),'MinPeakHeight',(max_dBFS-pow_min));
    
    if ~isempty(locs)
        if length(locs)>1
            [B,I]=sort(pks_fund{i},'descend')
            f_peaks=f(locs(I)+f_inf-1);
            for l=1:(length(f_peaks)-1)
                f_dif(l)=abs(f_peaks(l+1)-f_peaks(1));
            end            
            peak_loc=find(f_dif>50);
            if ~isempty(peak_loc)                    
                if length(peak_loc)>1
                    locs_fund{i}(2)=locs(find(pks_fund{i}==max(B(peak_loc+1))))+f_inf-1;
                    locs_fund{i}(1)=locs(I(1))+f_inf-1;
                else
                    locs_fund{i}(2)=locs(I(peak_loc+1))+f_inf-1;
                    locs_fund{i}(1)=locs(I(1))+f_inf-1;
                end
            else
                locs_fund{i}(1)=locs(I(1))+f_inf-1;
            end                  
        else
            locs_fund{i}=locs+f_inf-1;
        end
    end  
    clear f_dif
    i
end

% Sinusoides frecuencia fundamental
%[sinousoids3,frec_fund]=peakSta(locs_fund,dBFS,time_threshold,pow_threshold);
frec_fund=locs_fund;

for i=1:length(frec_fund)
    if length(frec_fund{i})>2
        frec_fund{i}=frec_fund{i}(1:2);
    end
end

%%% Vectores para plot tiempo
t=(1/fs)*linspace(0,length(s),size(frec_fund,2));


% Relación entre frecuencias
for i=1:length(frec_fund)
    f2{i}=f(frec_fund{i});
    if size(f2{i},2)==1
       rel_frec(1,i)=0; 
       f_x(i)=f2{i};%%Se le asigna una frecuencia x porque todavía no se sabe si es macho o hembra
       t_x(i)=1;
    elseif size(f2{i},2)==2       
       f_male(i)=max(f2{i});
       f_female(i)=min(f2{i});   
       rel_frec(1,i)=f_male(i)/f_female(i);
    else
       disp('error, se encontró más de una frecuencia fundamental')        
    end
end

ind_m=find(f_male~=0);
ind_f=find(f_female~=0);
f_male_mean=mean(f_male(ind_m));
f_female_mean=mean(f_female(ind_f));

for i=1:length(f_x) %%A cada frecuencia incognito se la asigna un genero mediante la resta con el promedio
    if t_x(i)==1
        if (abs(f_male_mean-f_x(i)))<(abs(f_female_mean-f_x(i)))
            f_male(i)=f_x(i);
        else
            f_female(i)=f_x(i);
        end
    end
end

% Gráfica frecuencia fundamental
figure; 
subplot(2,1,1); hold on
ind_m=find(f_male~=0);
ind_f=find(f_female~=0);
plot(t(ind_m),f_male(ind_m),'b');
plot(t(ind_f),f_female(ind_f),'r');
axis([0 t(end) f_i f_s])
hold off

% Gráfica relaciónes 
subplot(2,1,2); hold on
t_rel=linspace(0,t(end),length(rel_frec));
plot(t_rel,rel_frec,'k')
line=(1.5+1.5*rel_threshold)*ones(length(rel_frec),1);
plot(t_rel,line,'g')
line=(1.5-1.5*rel_threshold)*ones(length(rel_frec),1);
plot(t_rel,line,'g')
line=(1.33333+1.333333*rel_threshold)*ones(length(rel_frec),1);
plot(t_rel,line,'g')
line=(1.33333-1.333333*rel_threshold)*ones(length(rel_frec),1);
plot(t_rel,line,'g')
axis([0 t(end) 1.2 1.6]); hold off


f_conv=0;
cont=0;
j=1;
f_x=0;
t_x=0;
for i=1:length(rel_frec)
    if (rel_frec(i)<(1.5+1.5*rel_threshold) && rel_frec(i)>(1.5-1.5*rel_threshold)) || (rel_frec(i)<(1.33333+1.333333*rel_threshold) && rel_frec(i)>(1.333333-1.333333*rel_threshold))
       if f_conv==0
           f_conv=1;
           t_0(j)=t_rel(i);
           i_0(j)=i;
       else
           if f_x==1  
               t_x=0;
               f_x=0;
               cont=0;
           end
       end
    else
        if f_conv==1
           if (t(i)-t_0(j))<1
               f_conv=0;
               cont=0;
               f_x=0;
               t_0(j)=[];
               i_0(j)=[];
               clear t_x_0
               if isempty(t_0)
                   clear t_0 i_0 
               end                 
           else
               if f_x==0
                  t_x_0=t_rel(i);
                  f_x=1;
               else 
                   t_x=t_rel(i)-t_x_0;
                   cont=cont+1;
                   if t_x>1%1s
                      t_f(j)=t(i-cont);
                      i_f(j)=i-cont;
                      cont=0;
%                       clear t_x t_x_0
                      f_conv=0;
                      f_x=0;
                      j=j+1;              
                   end               
               end               
           end
        end
    end
    i
    if i==length(rel_frec) && f_conv==1
        t_f(j)=t_rel(i-cont);
        i_f(j)=i-cont;
    end
end


if exist('t_0','var') && exist('t_f','var')           
    for i=1:length(t_0)
        %%Gráfica limites convergencia
        subplot(2,1,1); hold on
        plot(t_0(i)*ones(800),linspace(0,800,800),'m');
        plot(t_f(i)*ones(800),linspace(0,800,800),'c');
        hold off;

        subplot(2,1,2); hold on
        plot(t_0(i)*ones(300),linspace(0,1.6,300),'m');
        plot(t_f(i)*ones(300),linspace(0,1.6,300),'c'); 
        hold off
        %%%Medidas de perturbación
        in_f(i)=ind_f(find(ind_f==i_0(i)));
        fin_f(i)=ind_f(find(ind_f==(i_f(i)-1)));
        in_m(i)=ind_m(find(ind_m==i_0(i)));
        fin_m(i)=ind_m(find(ind_m==(i_f(i)-1)));

        [ppq_mc(i),j_abs_m(i)]=perturbation(f_male(in_m(i):fin_m(i)),[],ppq_num);
        [ppq_fc(i),j_abs_f(i)]=perturbation(f_female(in_f(i):fin_f(i)),[],ppq_num);
        R_conv(i,:)=[ppq_mc(i)  ppq_fc(i) j_abs_m(i) j_abs_f(i) mean(f_male(i_0(i):i_f(i)-1)) mean(f_female(i_0(i):i_f(i)-1))]; 
        T_conv(i,:)=[(t_f(i)-t_0(i)) mean(rel_frec(i_0(i):i_f(i)))];
    end
else
    disp('no hubo convergencia') 
end

%Results
 
if exist('t_0','var') && exist('t_f','var') 
     if length(f_male(in_m:fin_m))>(length(ind_m)-20) || abs((ind_m(1)-in_m(1)))<(2*ppq_num)
        ppq_m=0;
        j_abs_m=0;
        Fm_m=0;
        Fm_mo=0;
    else
        [ppq_m,j_abs_m]=perturbation(f_male(ind_m((find(ind_m<=in_m(1))))),[],ppq_num);
        Fm_m=mean(f_male(ind_m((find(ind_m<=in_m(1))))));
        Fm_mo=mode(f_male(ind_m((find(ind_m<=in_m(1))))));;
    end

    if length(f_female(in_f:fin_f))>(length(ind_f)-20) || abs((ind_f(1)-in_f(1)))<(2*ppq_num)
        ppq_f=0;
        j_abs_f=0;
        Ff_m=0;
        Ff_mo=0;
    else
        [ppq_f,j_abs_f]=perturbation(f_female(ind_f(find(ind_f<=in_f(1)))),[],ppq_num);
        Ff_m=mean(f_female(ind_f(find(ind_f<=in_f(1)))));
        Ff_mo=mode(f_female(ind_f(find(ind_f<=in_f(1)))));
    end    
    %Prior to conv
    P=[ppq_m  ppq_f j_abs_m j_abs_f];      
    F=[Fm_m Ff_m Fm_mo Ff_mo];
    R0=[P F];        
    %Convergence
    R_conv=R_conv(1,:);
    T_conv=T_conv(1,:); 
    M_conv=[mode(f_male(in_m(1):fin_m(1)))  mode(f_female(in_f(1):fin_f(1)))];
    RC=[R_conv M_conv T_conv];
    %Results
    R=[R0 RC];    
else
    [ppq_m,j_abs_m]=perturbation(f_male(ind_m),[],ppq_num);
    [ppq_f,j_abs_f]=perturbation(f_female(ind_f),[],ppq_num);

    P=[ppq_m  ppq_f j_abs_m j_abs_f];
    M=[mode(f_male(ind_m)) mode(f_female(ind_f))];
    R0=[P mean(f_male(ind_m)) mean(f_female(ind_f)) mode(f_male(ind_m)) mode(f_female(ind_f))];
    R=[R0 zeros(1,10)]
    disp('no hubo convergencia')     
end

FileName
% 
 

