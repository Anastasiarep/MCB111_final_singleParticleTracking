%%
%diffusion eq

t=1:100:1000;
D1=30;
D2=300;
a=0:0.2:1;
r_sq=1:10:1E9;

p=zeros(length(t),length(r_sq));
p95=zeros(length(t),4,length(a));

%p=1-(a.*exp(-r_sq./(4.*D1.*t))+(1-a).*exp(-r_sq./(4.*D2.*t)));

for i=1:length(t)
    time=t(i);
    for j=1:length(a)
        alpha=a(j);
        prob=1-(alpha.*exp(-r_sq./(4.*D1.*time))+...
            (1-alpha).*exp(-r_sq./(4.*D2.*time)));
        [m,r] = min(abs(prob-0.95));
        p95(i,1,j)=r; %index
        p95(i,2,j)=time; %time
        p95(i,3,j)=r_sq(r); %r_sq
        p95(i,4,j)=prob(r); %0.95
    end    
end

%%

C={[204/255 0 0],[204/255 204/255 0],[0 204/255 0],...
    [0 102/255 204/255],[102/255 0 204/255],[204/255 0 102/255]};
hold on
for ii=1:6
    plot(p95(:,2,ii),p95(:,3,ii), 'Color', C{ii}, 'marker', '.')
end
set(gca, 'FontName', 'Cambria Math')
legend({'a=0';'a=0.2';'a=0.3';'a=0.4';'a=0.6';'a=1'},...
    'FontName', 'Cambria Math')
xlabel('Time', 'FontName', 'Cambria Math')
ylabel('r^{2} within which particle lives with p=0.95',...
    'FontName', 'Cambria Math')
title('0.95 probability r^{2} in time for different values of \alpha',...
    'FontName', 'Cambria Math')

x_saveas='rsqVtValpha.png';
cd 'C:\Users\anr612\Documents\1. Harvard18-\1. Courses\1. G1Spring-MCB111-Math\MCB111_FinalProject'
set(gcf, 'Color', 'w');
saveas(gcf, x_saveas);
