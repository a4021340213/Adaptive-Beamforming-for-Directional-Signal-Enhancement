function reconstructed_signal = noise_filter(noisy_signal,t)
m=0.5*numel(noisy_signal)+1;
n=0.5*numel(noisy_signal);
N=numel(t);
%% step 1 construct hankel matrix
for i=1:m
   
    for j=1:n
        H2(i,j)=noisy_signal(i+j-1);  %hankel matrix
    end
end

%% step 2 SVD
[U, sigma, V] = svd(H2);


%% step 3 correct mu value
total=sum(sum(sigma).^2);
sigma2=zeros(n,n);

for i=1:n
    
    for j=1:n
        if i==j
        sigma2(i,j)=sigma(i,j)^2+sqrt(total);
        end
    end
end
mu=sum(sigma2);

%% step 4 determine index k

AIC=zeros(1,n-1);
for d=1:n-1
        prod_term = prod(mu(d+1:n).^(1/(n-d)));
        mean_term = mean(mu(d+1:n));
        L_d = prod_term / mean_term;
        
        % Compute AIC(d)
        AIC_values(d) = -2 * N * (n - d) * log10(L_d) + 2 * d * (2 * n - d);
end

[~,k]=min(AIC_values);


%% step 5 ISVD

H_hat=zeros(size(H2));
for i=1:k
    H_hat=H_hat+U(:,i)*sigma(i,i)*V(:,i).';
end

%% step 6 Averaging method and getting result
reconstructed_signal=zeros(size(noisy_signal));
for i = 1:N
        % Calculate bounds l and h
        l = max(1, i - n + 1);
        h = min(n, i);
        
        % Compute the average as per the formula
        sum_term = 0;
        for j = l:h
            sum_term = sum_term + H_hat(i - j + 1, j);
        end
        
        % Normalize by the range
        reconstructed_signal(i) = sum_term / (h - l + 1);
end
    
 reconstructed_signal(N)=reconstructed_signal(N-1)

end