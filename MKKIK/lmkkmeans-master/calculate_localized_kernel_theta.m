function K_Theta = calculate_localized_kernel_theta(K, Theta)
K_Theta = zeros(size(K(:, :, 1)));
for m = 1:size(K, 3)
    K_Theta = K_Theta + (Theta(:, m) * Theta(:, m)') .* K(:, :, m);
end