%Smoothes a vector through a moving average
function [moving_avg] =movingaverage(data,DAYS_SEMI_WINDOW);
	DAYS=size(data,1);
	for n=1:DAYS
		if(n<DAYS_SEMI_WINDOW+1)
		   day_Window=[1:n+DAYS_SEMI_WINDOW]';
		elseif(n>DAYS-DAYS_SEMI_WINDOW)
		   day_Window=[n-DAYS_SEMI_WINDOW:DAYS]';
		else
			   day_Window=[n-DAYS_SEMI_WINDOW:n+DAYS_SEMI_WINDOW]';
		end

		moving_avg(n)=nanmean(data(day_Window));

	end
end
