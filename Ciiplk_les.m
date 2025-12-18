
% *****************************************************************************
%  RAÈUN KONSTITUCIJSKIH KOLIÈIN IN ÈLENOV TANGENTNE MATRIKE LESENEGA PREREZA
% *****************************************************************************
% uporabljen je elastièni zakon tako v nategu kot tlaku
% prirastek geometrijske deformacije je sestavljen in prirastka elastiène,
% temperaturne deformacij, deformacij lezenja in krèenja ter iz prirastka
% mehano-sorptivne deformacije


function [c11,c12,c22,defc_pr,Nc,Mc,D_w,sigma_w,Ds_w,Dc_w,Dms_w,alfa_w,D_wp,qu_w]=...\n
   	Ciiplk_les(koord_x,xg1,wg1,deps,dkappa,st_lamel,st_prrazd,b_lamele,zT_lamele_sp,...\n
      zT_lamele_zg,z,Dw,sigmaw,EtT,Et_T,EcT,Ec_T,temp_w,dmoist_w,dDt_w,krcenje_w,...\n
      lezenje_w,mod_cw,mehso_w,Dsw,Dcw,Dmsw,ft_T,fc_T,element,inkrement,iter,fid,cas,alfaw,Dwp,quw);



    global multip;
	c11=0;
	c12=0;
	c22=0;
	Nc=0;
	Mc=0;
	D_w=Dw;
    sigma_w=sigmaw;
	Ds_w=Dsw; %shrinkage
  	Dc_w=Dcw; %creep
   	Dms_w=Dmsw; %mehanosorptive
    D_wp=Dwp;
    qu_w=quw;
    alfa_w=alfaw;

   	% velikost èasovnega intervala (v dnevih)
    if inkrement > 1
        dcas=cas(inkrement)-cas(inkrement-1);
	else
        dcas=0;
	end

   
    Nc_w=0;
	Mc_w=0;
	c11_w=0;
	c12_w=0;
    c22_w=0;

    % sprememba elastiènega modula
	dEc=Ec_T-EcT;
    dEt=Et_T-EtT;

	% elastièna deformacija iz prejšnjega inkrementa
    Dwe=Dw-Dwp;

    % -------------------------------------------------
	%  Raèun napetosti v integracijskih toèkah prereza
	% -------------------------------------------------
    for i=1:sum(st_lamel)
        for p=1:st_prrazd(i)
         
            dsigma_dDw=zeros(size(xg1,2),size(xg1,2));

   	        for m=1:size(xg1,2)
                for n=1:size(xg1,2)
                    
                    
                    % raèun napetosti v integracijski toèki ob upoštevanju
                    % krèenja lesa
                    if krcenje_w==1
                        dDs_w(m,n,i,p)=multip(1)*0.00625*dmoist_w(m,n,i,p);    
                        
                        % celotna deformacija krèenja
                        Ds_w(m,n,i,p)=Dsw(m,n,i,p)+dDs_w(m,n,i,p);
                    else % brez upoštevanja krèenja
                        dDs_w(m,n,i,p)=0;
                    end
                    
                    % raèun napetosti v integracijski toèki ob upoštevanju
                    % lezenja lesa
                    % spremenljivka èas nastopa v [dnevih] (pretvorba [sekunda]->[dan])
                    if lezenje_w==1
                        % model A
                        if mod_cw==1
                            if sigmaw(m,n,i,p)>=0 % natezna napetost na zaèetku èasovnega inkrementa
                                if inkrement==1
                                   dDc_w(m,n,i,p)=0;  %dDc_w(m,n,i,p)=0.9927e-4;                
                                else
                                    alf1_T=multip(4)*-0.2719e-4;
                                    alf2_T=multip(5)*1.975e-2;
%                                       alf1_T=-1.83e-3;
%                                     alf2_T=2.35e-2;
                                    dDc_w(m,n,i,p)=sigmaw(m,n,i,p)*(alf1_T)*(exp(-alf2_T*cas(inkrement)/86400)-exp(-alf2_T*cas(inkrement-1)/86400));  
%                                           dDc_w(m,n,i,p)=sigmaw(m,n,i,p)*(-1.83e-3)*(exp(-2.35e-1*cas(inkrement)/86400)-exp(-2.35e-1*cas(inkrement-1)/86400));  
%                                     dDc_w(m,n,i,p)=sigmaw(m,n,i,p)*(-1.83e-3)*(exp(-2.35e-1*cas(inkrement)/86400)-exp(-2.35e-1*cas(inkrement-1)/86400));  
%                                    dDc_w(m,n,i,p)=sigmaw(m,n,i,p)*(-1.42788e-4)*(exp(-2.02961e-1*cas(inkrement)/86400)-exp(-2.02961e-1*cas(inkrement-1)/86400));  
%                                   dDc_w(m,n,i,p)=sigmaw(m,n,i,p)*(-1.83e-6)*(exp(-2.35e-5*cas(inkrement)/86400)-exp(-2.35e-5*cas(inkrement-1)/86400));
                                end    
                            else % tlaèna napetost v integracijski toèki (ustrezno razmerje kot v clanku)
                                if inkrement==1
                                 dDc_w(m,n,i,p)=0;    %dDc_w(m,n,i,p)=1.0080e-4;    
                                else    
                                    alf1_C=multip(2)*-0.1850e-4;
                                    alf2_C=multip(3)*2.017e-2;
%                                     alf1_C=-1.225e-3;
%                                     alf2_C=2.4e-2;
                                    dDc_w(m,n,i,p)=sigmaw(m,n,i,p)*(alf1_C)*(exp(-alf2_C*cas(inkrement)/86400)-exp(-alf2_C*cas(inkrement-1)/86400)); 
%                                     dDc_w(m,n,i,p)=sigmaw(m,n,i,p)*(-1.225e-3)*(exp(-2.399974684e-1*cas(inkrement)/86400)-exp(-2.399974684e-1*cas(inkrement-1)/86400));  
%                                     dDc_w(m,n,i,p)=sigmaw(m,n,i,p)*(-1.224935638e-3)*(exp(-2.399974684e-1*cas(inkrement)/86400)-exp(-2.399974684e-1*cas(inkrement-1)/86400));  
%                                     dDc_w(m,n,i,p)=sigmaw(m,n,i,p)*(-0.95577e-4)*(exp(-2.07277e-1*cas(inkrement)/86400)-exp(-2.0e-1*cas(inkrement-1)/86400));  
%                                    dDc_w(m,n,i,p)=sigmaw(m,n,i,p)*(-1.225e-6)*(exp(-2.4e-5*cas(inkrement)/86400)-exp(-2.4e-5*cas(inkrement-1)/86400));
                                end
                            end
                        end  %mod_cw==1 
 
                        %model B
                        if mod_cw==2
                            dDc_w(m,n,i,p)=sigmaw(m,n,i,p)*0.000018*...\n
                                sum(([0.0686 -0.0056 0.0716 0.0404 0.2073 0.5503]).*(1-exp(-(dcas/86400)./[0.01 0.1 1 10 100 5000])));
                        end  %mod_cw==2 
                                                    
                        % celotna deformacija lezenja
                        Dc_w(m,n,i,p)=Dcw(m,n,i,p)+dDc_w(m,n,i,p);
                        
                    else % brez upoštevanja lezenja
                        dDc_w(m,n,i,p)=0;
                        Dc_w(m,n,i,p)=Dcw(m,n,i,p)+dDc_w(m,n,i,p);
                    end
                    
                    % raèun napetosti v integracijski toèki ob upoštevanju
                    % mehano-sorptivne deformacije lesa
                    if mehso_w==1
                        if dmoist_w(m,n,i,p) >= 0 % pozitivna sprememba vlage v èasovnem inkrementu
                            dDms_w(m,n,i,p)=multip(6)*sigmaw(m,n,i,p)*0.00003*(1-exp(-1.6*abs(dmoist_w(m,n,i,p))));
                        else % vlaga se zmanjša
                            dDms_w(m,n,i,p)=multip(6)*sigmaw(m,n,i,p)*0.00003*(1-exp(-2.4*abs(dmoist_w(m,n,i,p))));
                        end    
                        % celotna mehano-sorptivna deformacija
                        Dms_w(m,n,i,p)=Dmsw(m,n,i,p)+dDms_w(m,n,i,p);
                        
                    else % brez upoštevanja mehano-sorptivne deformacije lesa
                        dDms_w(m,n,i,p)=0;
                    end
                  
                    % prirastek mehanske deformacije v integracijski toèki lesenega prereza
                    dDw=deps+z(m,n,i)*dkappa-dDt_w(m,n,i,p)-dDs_w(m,n,i,p)-dDc_w(m,n,i,p)-dDms_w(m,n,i,p);
		
	                % celotna mehanska deformacija v integracijski toèki 
                    D_w(m,n,i,p)=Dw(m,n,i,p)+dDw;

                   
                    
                    % tlaèno obmoèje konstitucijskega zakona lesa
                    if D_w(m,n,i,p) < 0 
                         % pomožna napetost na mestu integracijske toèke - sigma_trial
                         sigma_str=sigmaw(m,n,i,p)+dDw*Ec_T(m,n,i,p)+dEc(m,n,i,p)*Dwe(m,n,i,p);
                         % vpeljemo "relative stress"
                         ksi_str=sigma_str-quw(m,n,i,p);

                          % pomožna funkcija "f_trial"
                         f_trial=abs(ksi_str)-fc_T(m,n,i,p);

                        if f_trial <= 0 % elastièni korak
                              D_wp(m,n,i,p)=Dwp(m,n,i,p); % plastièna deformacija
                              dsigma_dDw(m,n)=Ec_T(m,n,i,p); % tangentni modul
                              sigma_w(m,n,i,p)=sigma_str; % napetost
                              qu_w(m,n,i,p)=quw(m,n,i,p); % zaostala napetost
                              alfa_w(m,n,i,p)=alfaw(m,n,i,p);
                               
                         else % (f_trial > 0) plastièni korak
                              % prirastek plastiènega konsistenènega parametra
                              Ep_T=15;
                              d_gama=f_trial/(Ec_T(m,n,i,p)+...\n
                                 Ec_T(m,n,i,p)*Ep_T/(Ec_T(m,n,i,p)-Ep_T));
                              alfa_w(m,n,i,p)=alfaw(m,n,i,p)+d_gama;
                              % plastièna deformacija                      
                              D_wp(m,n,i,p)=Dwp(m,n,i,p)+d_gama*sign(ksi_str);
                              % tangentni modul
                              dsigma_dDw(m,n)=Ep_T;
                              % napetost
                              sigma_w(m,n,i,p)=sigma_str-d_gama*Ec_T(m,n,i,p)*sign(ksi_str);                        
                              % zaostala napetost
                              qu_w(m,n,i,p)=quw(m,n,i,p)+...\n
                                d_gama*(Ec_T(m,n,i,p)*Ep_T/(Ec_T(m,n,i,p)-Ep_T))*sign(ksi_str); 
                        end % (f_trial > 0)  
                      
                    end
                    % natezno obmoèje konstitucijskega zakona lesa
                    if D_w(m,n,i,p) >= 0 
%                         sigma_w(m,n,i,p)=Et_T(m,n,i,p)*D_w(m,n,i,p);
% 				        dsigma_dDw(m,n)=Et_T(m,n,i,p);
                        
                        % pomožna napetost na mestu integracijske toèke - sigma_trial
                        sigma_str=sigmaw(m,n,i,p)+dDw*Et_T(m,n,i,p)+dEt(m,n,i,p)*Dwe(m,n,i,p);
                        % vpeljemo "relative stress"
                        ksi_str=sigma_str-quw(m,n,i,p);

                        % pomožna funkcija "f_trial"
                        f_trial=abs(ksi_str)-ft_T(m,n,i,p);
                        
                        if f_trial <= 0 % elastièni korak
                              D_wp(m,n,i,p)=Dwp(m,n,i,p); % plastièna deformacija
                              dsigma_dDw(m,n)=Et_T(m,n,i,p); % tangentni modul
                              sigma_w(m,n,i,p)=sigma_str; % napetost
                              qu_w(m,n,i,p)=quw(m,n,i,p); % zaostala napetost
                              alfa_w(m,n,i,p)=alfaw(m,n,i,p);
                               
                         else % (f_trial > 0) plastièni korak
                              % prirastek plastiènega konsistenènega parametra
                              Ep_T=15;
                              d_gama=f_trial/(Et_T(m,n,i,p)+...\n
                                Et_T(m,n,i,p)*Ep_T/(Et_T(m,n,i,p)-Ep_T));
                              alfa_w(m,n,i,p)=alfaw(m,n,i,p)+d_gama;
                              % plastièna deformacija                      
                              D_wp(m,n,i,p)=Dwp(m,n,i,p)+d_gama*sign(ksi_str);
                              % tangentni modul
                              dsigma_dDw(m,n)=Ep_T;
                              % napetost
                              sigma_w(m,n,i,p)=sigma_str-d_gama*Et_T(m,n,i,p)*sign(ksi_str);                        
                              % zaostala napetost
                              qu_w(m,n,i,p)=quw(m,n,i,p)+...\n
                                d_gama*(Et_T(m,n,i,p)*Ep_T/(Et_T(m,n,i,p)-Ep_T))*sign(ksi_str); 
                        end % (f_trial > 0) 
                        
                    end
            end %n=1:size(xg1,2)
         end %m=1:size(xg1,2)               
               
               
        % konstitucijska osna sila in upogibni moment obravnavanega lesenega preènega prereza
        Nc_w=Nc_w+b_lamele(i)/(4*st_prrazd(i))*(zT_lamele_sp(i)-zT_lamele_zg(i))/2*...\n
            sum(sum((wg1'*wg1).*(sigma_w(:,:,i,p))));
        Mc_w=Mc_w+b_lamele(i)/(4*st_prrazd(i))*(zT_lamele_sp(i)-zT_lamele_zg(i))/2*...\n
            sum(sum((wg1'*wg1).*(sigma_w(:,:,i,p)).*z(:,:,i)));
         
         
        % raèun èlenov konstitucijske togostne matrike lesenega preènega prereza
        c11_w=c11_w+b_lamele(i)/(4*st_prrazd(i))*(zT_lamele_sp(i)-zT_lamele_zg(i))/2*...\n
            sum(sum((wg1'*wg1).*dsigma_dDw));
        c12_w=c12_w+b_lamele(i)/(4*st_prrazd(i))*(zT_lamele_sp(i)-zT_lamele_zg(i))/2*...\n
            sum(sum((wg1'*wg1).*dsigma_dDw.*z(:,:,i)));
        c22_w=c22_w+b_lamele(i)/(4*st_prrazd(i))*(zT_lamele_sp(i)-zT_lamele_zg(i))/2*...\n
            sum(sum((wg1'*wg1).*dsigma_dDw.*z(:,:,i).^2));
        clear dsigma_dDw
	      
    end %p=1:st_prrazd(i)
end % i=1:sum(st_lamel)               
               
     
           

% ker se obravnava le polovica preènega prereza, moramo vsoto pomnožiti z 2
c11=2*c11_w;
c12=2*c12_w;
c22=2*c22_w;
% raèun konstitucijskih kolièin celotnega prereza 
Nc=2*Nc_w;
Mc=2*Mc_w;


%sestavljanje materialne tangentne matrike
C_l=[c11 c12
   c12 c22];
 
   
% definitnost tangentne konstitucijske matrike preènega prereza
if det(C_l) > 0 %& c11 > 0   
	% pozitivna definitnost
   defc_pr=1;
else
   % negativna definitnost
   defc_pr=-1;
end 
