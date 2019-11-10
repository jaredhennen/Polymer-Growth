function [lengths] = runnooverlappolymergrowth(int,na,nb,alpha,beta)
%Function to repeatedly run simulationg of non-overlapping polymer growth
%with monomers a and b with populations Na and Nb, angles alpha and beta,
%respectively.

lengths=zeros(int,1); % initialize lengths array
if (na + nb) > 50 %Check for large segment length to perform parallel for loop with - overhead of initializing not worthwhile with short segments
    parfor i=1:int
        lengths(i)=nooverlapstuckpolymergrowth(na,nb,alpha,beta);
    end
else
    for i=1:int
        lengths(i)=nooverlapstuckpolymergrowth(na,nb,alpha,beta);
    end
end

end