function [y w B] = RAIS_LCMV(matX,theta_i_hat,theta_s_hat,N,a,L,d)

a_s = exp(j*pi*sind(theta_s_hat).*d.');
a_i = exp(j*pi*sind(theta_i_hat).*d.');
R=eye(N);
w=zeros(N,1);
beta=0.998;

%% step 1
g = [1; 10^-6];


%% step2 for loop
    for l=1:99

        C = [a_s(:,l), a_i(:,l)];
        P=(1-C*inv(C'*C)*C');
        G=C*inv(C'*C)*g;
        R=beta*R+(1-beta).*matX(:,l)*matX(:,l)';
        eigenvalues = eig(R);
        max_eigenvalue = max(eigenvalues);
        mu=(2/max_eigenvalue)-0.1;
        w=P.*(eye(N)-mu*R)*w+G;
    end

    for l=1:L

        y(l) = w'*matX(:,l);
            for  i= 1 : numel(a(1,:))
                B(i,l) = 20*log10(abs(w'*a(:,i)));
            end
    end
end