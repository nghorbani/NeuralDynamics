function rate = find_rate(values, threshold, T)
set=0;rate=0;
for i = 1:length(values)
    if set == 0
        if values(i)>threshold
            set = 1;
            rate = rate + 1;
        end
    else
       if values(i)<threshold
           set = 0;
       end
    end
end
rate = rate / T;