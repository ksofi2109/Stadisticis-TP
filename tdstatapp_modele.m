% Travaux dirigés de Statistique Appliquée sous Matlab
% Script de départ
clear all       % Efface les variables en mémoire
clc             % Efface le contenu de la fenêtre de commande
format compact  % Supprime les sauts de ligne dans la fenêtre de commande
cas = 'B';      % Exercice choisi (cf. le switch sur sa valeur)
ALPHA = 0.05;   % Risque d'erreur de première espèce
chemin = '';    % Chemin pour les fichiers de données, par exemple : '' en local, '/home/esa/irivals/MATLAB2A/' en salle info

switch cas
    case 'A'
        display('Distributions de tailles (Galton)');
        noms = {'pères', 'mères', 'couples'};
        nom_fichier = sprintf('%sdataA.txt', chemin);
        x = load(nom_fichier);
        x(:, 3) = (x(:, 1) + 1.08*x(:,2))/2;
        
        for i = 1:3
            fprintf('\n*** %s ***\n', noms{i});
            tailles = x(:,i);
            
            % Tracer l'histogramme
            subplot(3,3,i)
            [eff,xb] = hist(tailles);
            plot(xb,eff)
            bar(xb,eff)
            
            % Tracer l'histogramme cumulé, etc.
            effcumu = cumsum(eff);
            subplot(3,3,i+3)
            bar(xb,effcumu)
            
            %Calcule des estimations biases
            esperance = mean(tailles);
            variance = var(tailles);        %Variance empirique S^2
            ecart_type = sqrt(variance);
            
            %Calcule des estimations no-biases
            variance_nb = var(tailles,1);
            
            %Coefficients
            coeff_d_asymetrie = skewness(tailles,0)       %0 si symetrie,3 si gaussien
            coeff_d_apltissement = kurtosis(tailles,0)
            
            %Testter le caractère gaussien (avec X^2)
            subplot(3,3,i+6)
            [h, p, stats] = chi2gof(tailles);
            qqplot(tailles);
        end
        
    otherwise
        display('PCRq sur des données de trisomie 21');
        nom_fichier = sprintf('%sdataB.txt', chemin);
        x = load(nom_fichier);

        % Moyenne des triplicats (i.e. moyennes des colonnes)
        HPRTmoy = mean(x(:,1:3),2);
        TRIOBPmoy = mean(x(:,4:6),2);
        MYLIPmoy  = mean(x(:,7:9),2);
        
        % Normaliser par rapport à HPRT (soustraction des ct)
        MYLIPnorm = MYLIPmoy-HPRTmoy;
        TRIOBPnorm = TRIOBPmoy-HPRTmoy;
        data = [TRIOBPnorm MYLIPnorm];
        
        % Tracer les boîtes à moustaches
        groupe = [zeros(10,1); ones(10,1)];
        figure(1)
        subplot(1,2,1)
        boxplot(TRIOBPnorm, groupe, 'labels', {'T21+', 'T21-'})
        title('TRIOP normé')
        
        subplot(1, 2, 2)
        boxplot(MYLIPnorm, groupe, 'labels', {'T21+', 'T21-'})
        title('MYLOP normé')
        
        % Analyse différentielle (comparaison des TS21+ et TS21-)
        noms = {'TRIOBP', 'MYLIP'};
        for i=1:2
            fprintf('\n*** Gène %s ***\n', noms{i});
            gene = data(:,i);
            geneplus = gene(1:10);
            genemoins = gene(11:20);
            
            %Test de normalité
            [hp, pp, wp] = swtest(geneplus);
            [hm, pm, wm] = swtest(genemoins);
            
            figure(2)
            subplot(2,2,i)
            qqplot(geneplus)
            titre = sprintf('%s Shapiro-Wilk T21+ : p= %.2g (W = %.2f)',noms{i}, pp, wp);
            title(titre)
            
            figure(2)
            subplot(2,2,i+2)
            qqplot(genemoins)
            titre = sprintf('%s Shapiro-Wilk T21- : p= %.2g (W = %.2f)',noms{i}, pm, wm);
            title(titre)
            
            if (hp == 0 & hm == 0)
                [hv, pv, civ, statsv] = vartest2(geneplus, genemoins);
                if hv == 0
                    [he, pe, cie, statse] = ttest2(geneplus, genemoins, 'vartype', 'equal');
                else
                    [he, pe, cie, statse] = ttest2(geneplus, genemoins, 'vartype', 'unequal');
                end
            else
                [pw, hw, statsw] = ranksum(geneplus, genemoins);
                nboot = 10000;
                mplus = bootstrp(nboot, @mean, geneplus);
                mmoins = bootstrp(nboot, @mean, genemoins);
                d = mplus-mmoins;
                q = quantile(d, [0.025 0.975]);
                
                figure(3), clf, set(3, 'position', [100 100 500 500])
                hist(d)
                
            end
        end
end

