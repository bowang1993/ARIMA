function z=avspec(y,n,slide,plfg)
%
%      z=avspec(y,n,slide,plfg)
%
% computes an average of periodograms using an FFT of length n
% moves by slide between each successive periodogram
% z is average over batches on linear scale
% plfg >0 to plot 
% may be adapted with a specific window; check the code to see
%
ny=length(y);
nbat=floor((ny-n)/slide)+1;
disp(['no. batches = ',num2str(nbat)]);
nb=1;
ne=nb+n-1;
yav=zeros(1,floor(n/2));
for i=1:nbat
%	yin=hanning(n) .* y(nb:ne);
	yin=y(nb:ne);
	yout=abs(fft(yin));
	yout=yout(1:floor(n/2));
	yav=yav+(yout.*yout)/(2*pi*n);
	nb=nb+slide;
	ne=ne+slide;
end
z=yav/nbat;
if plfg
    semilogy(log(yav/nbat));
    title(['Average of ',int2str(nbat),' periodograms of length ',int2str(n)]);
end
