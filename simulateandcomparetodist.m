function [] = simulateandcomparetodist(nit,na,nb,alpha,beta,binfact)
%Function to calculate a computational distribution for polymer
%growth which allows for overlapping, then compare it to the expected
%analytical distribution.
%Computational distribution determined from nit iterations
%Polymer consists of monomers a and b with populations Na and Nb, bond
%angles alpha and beta, respectively

a = (2*sin(alpha/2))^2*na; b = (2*sin(beta/2))^2*nb; %effective squared lengths of each component for analytical distribution
maxv = sqrt(a+b)*3; %determine the max value to plot to, 3x the expected RMS value

lengths = simulatemanypolymerswoverlap(nit,na,nb,alpha,beta); %determine the computational distribution
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

end