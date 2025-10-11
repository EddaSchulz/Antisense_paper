function calc_error = fit_k1_off(p, t_data, X_data, kdeg,Y0)
    p=exp(p);
    [~, Y] = ode45(@(t, Y) model_diff(t, Y, p, kdeg), t_data, Y0);
    X_model = Y;  % Extract protein values
    calc_error = sum((X_model - X_data).^2);
end
