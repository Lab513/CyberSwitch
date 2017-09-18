function cost = cost_fcn(exp_data, sim_data, duration, weight)

    %cost= weight*sum(((mean(exp_data,2) - sim_data) ./ (mean(exp_data,2)+ sim_data)).^2) / duration;
    cost= weight*sum((mean(exp_data,2) - sim_data).^2) / duration;
    
end

