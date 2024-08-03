function [t_all, x_all] = rk4(f, t_span, x0)

n_steps = numel(t_span);
n_states = numel(x0);
x_all = NaN(n_steps,n_states);
x_all(1,:) = x0;
t_all = reshape(t_span,n_steps,1);

for i = 1:n_steps - 1
    h = t_span(i+1) - t_span(i);
    x = x_all(i,:)';
    t = t_span(i);
    
    k1 = f(t, x);
    k2 = f(t + h/2, x + h*k1/2);
    k3 = f(t + h/2, x + h*k2/2);
    k4 = f(t + h, x + h*k3);
    
    x_all(i+1,:) = x + h*(k1 + 2*k2 + 2*k3 + k4)/6;
    
    
end

end