clear all
b=0.8;%i_rf
c=0:0.01:1; %i_dc
x=length(c);

h=6.626*10.^-34/(2*pi);%Planck's constant
f=(2*pi)*2*10.^9;%rf frequency in Hz
e=1.6*10.^-19;%electron charge
Ic=3.3*10.^-6;%critical current
R=30;%constant (background) resistance
tau=(2*e*R*Ic)/(h);

% Code for choosing i_dc for the odd steps
avgphidot=zeros(1,x);%average voltage
phi_initial=0;
phidot_initial=0;
    
for m=1:x
        
    step=0.0003*10.^-9;
    t=0:0.0003*10.^-9:6*10.^-9; %t in seconds
    k=length(t);
        
    phi=zeros(1,k);
    phi(1)=phi_initial;
    phidotnew=zeros(1,k);
    phidotnew(1)=phidot_initial;
      
    %RK4 method
    for i=1:k-1
            
        phidot=@(t,phi)(tau*(b*sin(f*t)-sin(phi)+c(m)));
            
        k1 = phidot(t(i),phi(i));
        k2 = phidot(t(i)+0.5*step,phi(i)+0.5*k1*step);
        k3 = phidot(t(i)+0.5*step,phi(i)+0.5*k2*step);
        k4 = phidot(t(i)+step,phi(i)+k3*step);

        phi(i+1) = phi(i)+((k1+2*k2+2*k3+k4)/6)*step;
        phidotnew(i+1)=phidot(t(i+1),phi(i+1));
            
    end
        
    avgphidot(1,m)=mean(phidotnew)/(f);
    phi_initial=phi(k);
    phidot_initial=phidotnew(k);
end
   
plot(avgphidot(1,1:x),c)
xlabel('{\it V/[hf/2e]}', 'interpreter', 'tex')
ylabel('{\it i_{dc}}', 'interpreter', 'tex')
set(gca,'TickDir','out');
xlim([0,7])