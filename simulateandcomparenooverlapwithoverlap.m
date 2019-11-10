function [] = simulateandcomparenooverlapwithoverlap(nit,na,nb,alpha,beta,binfact)
%Function to compare chain length distributions of overlap not allowed and
%overlap allowed cases of polymer growth.
%Computational distribution determined from nit iterations
%Polymer consists of monomers a and b with populations Na and Nb, bond
%angles alpha and beta, respectively

lengths = runnooverlappolymergrowth(nit,na,nb,alpha,beta);

a = (2*sin(alpha/2))^2*na; b = (2*sin(beta/2))^2*nb; %effective squared lengths of each component for analytical distribution
maxv = max(lengths); %determine the max value to plot to, using maximum of computational distribution

h=histogram(lengths,(0:binfact:maxv),'Normalization','Probability'); %plot a histogram of computational distribution

theorydist = @ (x,a,b) 2*x./(a+b).*exp(-x.^2/(a+b)); %define expected analytical distribution

hold on;
x=(0:maxv/1000:maxv); %initialize chain length array for analytical distribution plot
p=plot(x,theorydist(x,a,b)*binfact,'LineWidth',3); %overplot analytical distribution
%make plot look nice
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold');
xlabel('||\Gamma|| (bond lengths)','FontWeight','bold','FontSize',12);
ylabel('P(||\Gamma||)','FontWeight','bold','FontSize',12);
hold off;