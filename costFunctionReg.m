function [J, grad] = costFunctionReg(theta, X, y, lambda)
%COSTFUNCTIONREG Compute cost and gradient for logistic regression with regularization
%   J = COSTFUNCTIONREG(theta, X, y, lambda) computes the cost of using
%   theta as the parameter for regularized logistic regression and the
%   gradient of the cost w.r.t. to the parameters. 

% Initialize some useful values
m = length(y); % number of training examples

% You need to return the following variables correctly 
J = 0;
grad = zeros(size(theta));

% Compute the cost of a particular choice of theta.
% You should set J to the cost.
% Compute the partial derivatives and set grad to the partial
% derivatives of the cost w.r.t. each parameter in theta


for i = 1:m
    h_theta_X(i) = sigmoid((theta'*X(i,:)'));
    J(i) = (-y(i)*log(h_theta_X(i)))-((1-y(i))*log(1-h_theta_X(i)));
end

for j = 2:size(theta)
    thetaJ(j) = theta(j)^2;
end

J = ((1/(m))*sum(J))+((lambda/(2*m))*sum(thetaJ));



for i = 1:m
    h_theta_X(i) = sigmoid((theta'*X(i,:)'));
    grad_temp(1,i) = (h_theta_X(i)-y(i))*X(i,1); %%% Supply a temporary variable for the gradient calculation for each hypothesis value
end

grad(1) = (1/m)*sum(grad_temp(1,:));

for j = 2:size(theta,1)
    for i = 1:m
        h_theta_X(i) = sigmoid((theta'*X(i,:)'));
        grad_temp(j,i) = (h_theta_X(i)-y(i))*X(i,j);
    end
    grad(j) = ((1/m)*(sum(grad_temp(j,:))))+((lambda/m)*theta(j));
end




% =============================================================

end
