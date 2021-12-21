function [label,pmax] = clasificar_one_bayes(x,betas,modelos,nclases,priors_groups)
logZ = -1e10;
paux = logZ*ones(1,nclases);
for i=1:length(betas), %loop sobre grupos
    W = betas{i};
    model = modelos{i};
    p_aux = logZ;
    for j=1:nclases % loop sobre clases
        if model.probs(j)~=0,
            p_aux = logsum(p_aux, log(model.probs(j)) + log(mvnpdf(x*W,model.means(:,j)',squeeze(model.sigmas(j,:,:)))));
        else
            p_aux = logsum(p_aux,logZ);
        end
    end
    for j=1:nclases % loop sobre clases
        if model.probs(j)~=0,
            paux(j) = logsum(paux(j), log(model.probs(j)) - p_aux + log(mvnpdf(x*W,model.means(:,j)',squeeze(model.sigmas(j,:,:)))));
        else
            paux(j) = logsum(paux(j),logZ);
        end
    end

end
label = argmax(paux);
pmax = exp(paux(label));