function x_B = RTN_to_ECI(x_rel, x_ref)

i_hat = x_ref(1:3)/norm(x_ref(1:3));
k_hat = cross(x_ref(1:3),x_ref(4:6))/norm(cross(x_ref(1:3),x_ref(4:6)));
j_hat = cross(k_hat, i_hat);

T = [i_hat'; j_hat'; k_hat']; % ECI to LVLH transform
r_rel = T'*x_rel(1:3); % Convert back to ECI frame
v_rel = T'*x_rel(4:6);
x_rel = [r_rel; v_rel];

r_A = x_ref(1:3);
v_A = x_ref(4:6);

r_B = x_rel(1:3) + r_A;

omega = cross(r_A, v_A)/(norm(r_A)^2);
v_B = x_rel(4:6) + v_A + cross(omega,x_rel(1:3));

x_B = [r_B; v_B];

end