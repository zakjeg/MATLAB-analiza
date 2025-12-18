
% *********************************************************************
%  NELINEARNA ANALIZA LESENIH KONSTRUKCIJ V SPREMENLJIVEM OKOLJU Z MKE
% *********************************************************************
% vgrajena bisekcija za iskanje limitnih stanj

% branje podatkov o geometriji konstrukcije iz datoteke PODATKI.M
Podatki

% zaèetni inkrement je enak 1
if zacetek == 0

clear all

disp(sprintf('\n *********************************************************************'))
disp(sprintf('          FAKULTETA ZA GRADBENIŠTVO IN GEODEZIJO, Ljubljana'))
disp(sprintf('              Katedra za masivne in lesene konstrukcije'))
disp(sprintf('\n     ---------- Leseni nosilci v spremenljivem okolju ----------'))
disp(sprintf('                      Izdelal: Sebastjan Bratina'))
disp(sprintf(' *********************************************************************\n'))



% branje podatkov o geometriji konstrukcije iz datoteke PODATKI.M
Podatki

%število vozlišè konstrukcije
st_vozlisc=size(voz,1);

%število elementov konstrukcije
st_elementov=size(pq,1);

% število razliènih preènih prerezov
if isempty(pr)
   st_prerezov=1;
else   
	st_prerezov=max(pr(:,2));
end
% število razliènih temperaturnih polj
if isempty(tpr)
   st_tpolj=1;
else   
	st_tpolj=max(tpr(:,2));
end


fid = fopen('izpis.out','wt');

%zapis na datoteko IZPIS
fprintf(fid,' *********************************************************************');
fprintf(fid,'\n          FAKULTETA ZA GRADBENISTVO IN GEODEZIJO, Ljubljana');
fprintf(fid,'\n              Katedra za masivne in lesene konstrukcije');
fprintf(fid,'\n\n     ---------- Leseni nosilci v spremenljivem okolju ----------');
fprintf(fid,'\n                     Izdelal: Sebastjan Bratina');
fprintf(fid,'\n *********************************************************************\n');
fprintf(fid,' -> Stevilo vozlisc konstrukcije: %i\n',st_vozlisc);
fprintf(fid,' -> Stevilo elementov konstrukcije: %i\n',st_elementov);
fprintf(fid,' -> Stevilo razlicnih precnih prerezov: %i\n',st_prerezov);
fprintf(fid,' -> Stevilo razlicnih temperaturnih polj: %i\n',st_tpolj);


fprintf(fid,'\n  Upoštevani deformacijski prispevki:');
if krcenje_w==1
   fprintf(fid,'\n -> krèenje lesa');
end
if lezenje_w==1
   if mod_cw==1
       fprintf(fid,'\n -> lezenje lesa, model A');
   else
       fprintf(fid,'\n -> lezenje lesa, model b');
   end  
end
if mehso_w==1
   fprintf(fid,'\n -> mehano-sorptivna def.\n');
end  



% podpore na konstrukciji
%disp(sprintf('\n\n\n *********************** PODPORE KONSTRUKCIJE ***********************'))
%disp(sprintf('\n%10s %11s %15s %15s','vozlišèe','smer X','smer Z','smer Y'))
%disp(' ====================================================================')
izpis_pod=zeros(1,3*st_vozlisc);
for i=1:st_vozlisc
   for j=((i-1)*3+1):3*i
      for k=1:size(fix,2)
         if fix(k)==j
            izpis_pod(j)=1;
         end   
      end
   end 
   izp(i,:)=izpis_pod(((i-1)*3+1):3*i);
end


disp(sprintf('\n\n -> Število vozlišè konstrukcije: %i',st_vozlisc))
disp(sprintf(' -> Število elementov konstrukcije: %i',st_elementov))
disp(sprintf(' -> Število razliènih preènih prerezov: %i',st_prerezov))

% najveèje število sprošèenih prostostnih stopenj elementa
for element=1:st_elementov
	if isempty(spr)
	else
		for i=1:size(spr,1)
			if spr(i,1)==element
	   		s_p(i)=size(find(spr(i,2:size(spr,2))==0),2);
	      end
   	end   
   end   
end   



	% branje podatkov o nespremenljivi obtežbi konstrukcije iz datotek OBTEZBA_NS.M
   obtezba_ns
   % branje podatkov o spremenljivi obtežbi konstrukcije iz datotek OBTEZBA_SP.M
   obtezba_sp

	% matrika nespremenljive in spremenljive vozlišène obtežbe konstrukcije v GKS  
   obtn_voz=zeros(st_vozlisc,3);
   obtsp_voz=zeros(st_vozlisc,3);
	for vozlisce=1:st_vozlisc
	   for i=1:size(obtn_X_voz,1)
	      if obtn_X_voz(i,1)==vozlisce
		      obtn_voz(vozlisce,1)=obtn_X_voz(i,2);
		   end 
      end 
      for i=1:size(obtsp_X_voz,1)
	      if obtsp_X_voz(i,1)==vozlisce
		      obtsp_voz(vozlisce,1)=obtsp_X_voz(i,2);
		   end 
      end
      
	   for i=1:size(obtn_Y_voz,1)
	      if obtn_Y_voz(i,1)==vozlisce
		      obtn_voz(vozlisce,3)=obtn_Y_voz(i,2);
		   end 
      end  
      for i=1:size(obtsp_Y_voz,1)
	      if obtsp_Y_voz(i,1)==vozlisce
		      obtsp_voz(vozlisce,3)=obtsp_Y_voz(i,2);
		   end 
      end
      
	   for i=1:size(obtn_Z_voz,1)
	      if obtn_Z_voz(i,1)==vozlisce
		      obtn_voz(vozlisce,2)=obtn_Z_voz(i,2);
		   end 
      end
      for i=1:size(obtsp_Z_voz,1)
	      if obtsp_Z_voz(i,1)==vozlisce
		      obtsp_voz(vozlisce,2)=obtsp_Z_voz(i,2);
		   end 
   	end
	end 
   
   % sprostitve elementov (zapis na datoteko)
   spr_zap=ones(st_elementov,6);
   
	% matrika nespremenljive in spremenljive obtežbe na elementih v LKS
   obtn_elem=zeros(st_elementov,6);
   obtsp_elem=zeros(st_elementov,6);
	for element=1:st_elementov
	  	for i=1:size(obtn_x_el,1)
	     	if obtn_x_el(i,1)==element
	         obtn_elem(element,1:2)=obtn_x_el(i,2:3);
	      end 
      end    
      for i=1:size(obtsp_x_el,1)
	     	if obtsp_x_el(i,1)==element
	         obtsp_elem(element,1:2)=obtsp_x_el(i,2:3);
	      end 
      end
      
	   for i=1:size(obtn_y_el,1)
	      if obtn_y_el(i,1)==element
	         obtn_elem(element,5:6)=obtn_y_el(i,2:3);
	      end
      end    
      for i=1:size(obtsp_y_el,1)
	      if obtsp_y_el(i,1)==element
	         obtsp_elem(element,5:6)=obtsp_y_el(i,2:3);
	      end
      end
      
	   for i=1:size(obtn_z_el,1)
	      if obtn_z_el(i,1)==element
	         obtn_elem(element,3:4)=obtn_z_el(i,2:3);
	      end
      end
      for i=1:size(obtsp_z_el,1)
	      if obtsp_z_el(i,1)==element
	         obtsp_elem(element,3:4)=obtsp_z_el(i,2:3);
	      end
      end
      
      
      % sprostitve elementov
      if isempty(spr)
	      sprost(element,:)=0;
      else
         sprost(element,max(s_p))=0;
			for i=1:size(spr,1)            
            if spr(i,1)==element
               k=1:6;
               k(find(spr(i,2:size(spr,2))))=[];
               for m=1:size(k,2)
                  sprost(element,m)=k(m);
                  spr_zap(element,sprost(element,m))=0;
               end   
            end   
	   	end   
      end  
   end 
   

  
   % izpis nespremenljive vozlišène obtežbe konstrukcije
	%disp(sprintf('\n\n ******* NESPREMENLJIVA VOZLIŠÈNA OBTEŽBA KONSTRUKCIJE (Kn,m) *******'))
	%disp(sprintf('\n%10s %11s %15s %15s','vozlišèe','sila X','sila Z','moment Y'))
	%disp(' ====================================================================')
	%for i=1:st_vozlisc
	%   disp(sprintf('%6i %14.1f %15.1f %15.1f',i,[obtn_voz(i,1) obtn_voz(i,2) obtn_voz(i,3)]'))
	%	disp(' --------------------------------------------------------------------')
   %end
   % izpis spremenljive vozlišène obtežbe konstrukcije
	%disp(sprintf('\n\n ******** SPREMENLJIVA VOZLIŠÈNA OBTEŽBA KONSTRUKCIJE (Kn,m) ********'))
	%disp(sprintf('\n%10s %11s %15s %15s','vozlišèe','sila X','sila Z','moment Y'))
	%disp(' ====================================================================')
	%for i=1:st_vozlisc
	%   disp(sprintf('%6i %14.1f %15.1f %15.1f',i,[obtsp_voz(i,1) obtsp_voz(i,2) obtsp_voz(i,3)]'))
	%	disp(' --------------------------------------------------------------------')
   %end
   
	% izpis nespremenljive obtežbe na elementih
	%disp(sprintf('\n\n ******* NESPREMENLJIVA OBTEŽBA NA ELEMENTIH KONSTRUKCIJE (Kn,m) *******'))
	%disp(sprintf('\n%26s %20s %20s','obtežba x','obtežba z','obtežba y'))
	%disp(sprintf('\n%9s %9s %9s %9s %9s %10s %9s','element','qx(0)','qx(L)','qz(0)','qz(L)','my(0)','my(L)'))
	%disp(' ========================================================================')
	%for i=1:st_elementov
	%   disp(sprintf('%6i %11.2f %9.2f %10.2f %9.2f %10.2f %9.2f',i,obtn_elem(i,:)'))
	%	disp(' ------------------------------------------------------------------------')
  % end
   % izpis spremenljive obtežbe na elementih
	%disp(sprintf('\n\n ******** SPREMENLJIVA OBTEŽBA NA ELEMENTIH KONSTRUKCIJE (Kn,m) ********'))
	%disp(sprintf('\n%26s %20s %20s','obtežba x','obtežba z','obtežba y'))
	%disp(sprintf('\n%9s %9s %9s %9s %9s %10s %9s','element','qx(0)','qx(L)','qz(0)','qz(L)','my(0)','my(L)'))
	%disp(' ========================================================================')
	%for i=1:st_elementov
	 %  disp(sprintf('%6i %11.2f %9.2f %10.2f %9.2f %10.2f %9.2f',i,obtsp_elem(i,:)'))
		%disp(' ------------------------------------------------------------------------')
	%end



% **************************************
%  RAÈUN KARAKTERISTIK PREÈNEGA PREREZA
% **************************************

% preberejo se podatki, ki so na datoteki PODATKI1.M
preberi
   
   
   
% merjenje raèunskega èasa   
tm=cputime;

for p=1:st_prerezov
   % raèun karakteristik lamel p-tega preènega prereza 
   % funkcija je zapisana v datoteki LAMELE1.M
	[ploscina(p),b_lam,zT_lamzg,zT_lamsp,zT(p),ztez_lok]=lamele1(st_podprerezov(p),...\n
      visina(p,:),b_sp(p,:),b_zg(p,:),st_lamel(p,:),p);
   if p==1
    	b_lamele=b_lam;
    	zT_lamele_zg=zT_lamzg;
      zT_lamele_sp=zT_lamsp;
      zT_lok=ztez_lok; 
    else
      b_lamele=([b_lamele,b_lam]);
   	zT_lamele_zg=([zT_lamele_zg,zT_lamzg]);
      zT_lamele_sp=([zT_lamele_sp,zT_lamsp]);
      zT_lok=([zT_lok,ztez_lok]); 
   end
   
   % Raèun vzrajnostnega momenta p_tega preènega prereza
   Iy(p)=Iy_prereza(st_podprerezov,p,b_sp(p,:),b_zg(p,:),visina(p,:),zT(p),zT_lok);
   
   k=1;
	% število preènih razdelitev posameznih podprerezov
	for i=1:size(st_lamel(p,:),2)
   	for j=1:st_lamel(p,i)
      	prrazd(p,k)=st_prrazd(p,i);
	      k=k+1;
   	end
	end
end

% število lamel posameznega, p-tega preènega prereza
p=1:st_prerezov;
if size(st_lamel,2)==1
   lamele=st_lamel;
else
   lamele(p)=sum(st_lamel(p,:)');
end




%raèun dolžin elementov (v centimetrih)
i=1:size(pq,1);
dolzina(i)=100*sqrt((voz(pq(i,3),2)-voz(pq(i,2),2)).^2+...\n
                     (voz(pq(i,3),3)-voz(pq(i,2),3)).^2);


% zapis karakteristik precnih prerezov na datoteko IZPIS
fprintf(fid,'\n\n                      ****************************\n');
fprintf(fid,'                       PODATKI O PRECNIH PREREZIH \n');  
fprintf(fid,'                      ****************************\n');

for p=1:st_prerezov
   % temperaturno polje obravnavanega preènega prereza
	if isempty(tpr)
   	polje=1;
	else   
		polje=tpr(p,2);
	end
	fprintf(fid,'\n\n   ******************** %i. PRECNI PREREZ (enote: kN,cm) ********************\n',p);
   fprintf(fid,'\n    zT = %10.2f\n',zT(p));
   fprintf(fid,'     A = %10.2f\n',ploscina(p));
   fprintf(fid,'   A_s = %10.2f\n',G_ploscina(p));
   fprintf(fid,'    Iy = %10.1f\n',Iy(p));
   fprintf(fid,'    temp.polje: %i.\n',polje); 
   fprintf(fid,' ========================================================================================\n\n');
  	fprintf(fid,'%13s %8s %9s %10s %10s %10s','LES:','Et','Ec','Eg','ft','fc');
   fprintf(fid,'\n               --------------------------------------------------------------------------\n');
   fprintf(fid,'%23.0f %9.0f %10.0f %10.2f %10.2f\n',[E_t0(p) E_c0(p) G_modul(p) f_t0(p) f_c0(p)]);
   fprintf(fid,'\n\n%12s %8s %9s %7s %17s %19s\n','PODPREREZ','b_zg','b_sp','h','stevilo lamel','stevilo vert.del.');
   fprintf(fid,' -----------------------------------------------------------------------------\n');
   fprintf(fid,'%8i %12.1f %9.1f %9.1f %9i %16i\n',[(1:st_podprerezov(p))' b_zg(p,1:st_podprerezov(p))' b_sp(p,1:st_podprerezov(p))' visina(p,1:st_podprerezov(p))' st_lamel(p,1:st_podprerezov(p))' prrazd(p,1:st_podprerezov(p))']');
   fprintf(fid,' =============================================================================\n\n');
   
end

% zapis karakteristik vozlišè na datoteko IZPIS
fprintf(fid,'\n\n\n                       ****************************\n');
fprintf(fid,'                        PODATKI O VOZLISCIH KONSTR. \n');  
fprintf(fid,'                       ****************************\n');
fprintf(fid,'\n%32s %44s\n','koordinate(m)','podpore (1-preprecen, 0-prosti premik)');
fprintf(fid,'%10s %9s %11s %16s %11s %11s\n','vozlisce','X','Z','smer X','smer Z','smer Y');
fprintf(fid,' ==============================================================================\n');
fprintf(fid,'%7.0f %14.2f %11.2f %12i %11i %10i\n',[voz izp]');
fprintf(fid,' ==============================================================================\n');


% zapis karakteristik elementov na datoteko IZPIS
%preèni prerez elementa
if isempty(pr)
   pre(1:st_elementov)=1;
else   
   pre(1:st_elementov)=pr(1:st_elementov,2);
end

fprintf(fid,'\n\n\n                      ****************************\n');
fprintf(fid,'                       PODATKI O ELEMENTIH KONSTR. \n');  
fprintf(fid,'                      ****************************\n');
fprintf(fid,'\n%21s %11s %12s %12s\n','zacetno','koncno','dolzina','precni');
fprintf(fid,'%10s %11s %11s %9s %14s %16s\n','element','vozlisce','vozlisce','(m)','prerez','sprostitve');
fprintf(fid,' ================================================================================\n');
fprintf(fid,'%7i %10i %11i %14.2f %11i %7i %2i %2i %2i %2i %2i\n', [pq';dolzina./100;pre;spr_zap']);
fprintf(fid,' ================================================================================\n');


% zapis nespremenljive vozlišène obtežbe konstrukcije na datoteko IZPIS
fprintf(fid,'\n\n\n   *******************************************************\n');
fprintf(fid,'    NESPREMENLJIVA VOZLISCNA OBTEZBA KONSTRUKCIJE (kN, m)\n');  
fprintf(fid,'   *******************************************************\n');
fprintf(fid,'\n%10s %11s %15s %16s\n','vozlisce','sila X','sila Z','moment Y');
fprintf(fid,' ==================================================================\n');
fprintf(fid,'%6i %14.1f %15.1f %15.1f\n',[(1:st_vozlisc)' obtn_voz]');
fprintf(fid,' ==================================================================\n');

% zapis spremenljive vozlišène obtežbe konstrukcije na datoteko IZPIS
fprintf(fid,'\n\n\n     *****************************************************\n');
fprintf(fid,'      SPREMENLJIVA VOZLISCNA OBTEZBA KONSTRUKCIJE (kN, m)\n');  
fprintf(fid,'     *****************************************************\n');
fprintf(fid,'\n%10s %11s %15s %16s\n','vozlisce','sila X','sila Z','moment Y');
fprintf(fid,' ==================================================================\n');
fprintf(fid,'%6i %14.1f %15.1f %15.1f\n',[(1:st_vozlisc)' obtsp_voz]');
fprintf(fid,' ==================================================================\n');
	
	
% zapis nespremenljive obtežbe na elementih konstrukcije na datoteko IZPIS
fprintf(fid,'\n\n\n         ****************************************************\n');
fprintf(fid,'          NESPREMENLJIVA OBTEZBA NA ELEMENTIH KOSTR. (kN, m)\n');  
fprintf(fid,'         ****************************************************\n');
fprintf(fid,'\n%26s %20s %20s\n','obtezba x','obtezba z','obtezba y');
fprintf(fid,'%9s %9s %9s %9s %9s %10s %9s\n','element','qx(0)','qx(L)','qz(0)','qz(L)','my(0)','my(L)');
fprintf(fid,' ========================================================================\n');
fprintf(fid,'%6i %11.2f %9.2f %10.2f %9.2f %10.2f %9.2f\n',[(1:st_elementov)' obtn_elem]');
fprintf(fid,' ========================================================================\n');

% zapis spremenljive obtežbe na elementih konstrukcije na datoteko IZPIS
fprintf(fid,'\n\n\n          **************************************************\n');
fprintf(fid,'           SPREMENLJIVA OBTEZBA NA ELEMENTIH KOSTR. (kN, m)\n');  
fprintf(fid,'          **************************************************\n');
fprintf(fid,'\n%26s %20s %20s\n','obtezba x','obtezba z','obtezba y');
fprintf(fid,'%9s %9s %9s %9s %9s %10s %9s\n','element','qx(0)','qx(L)','qz(0)','qz(L)','my(0)','my(L)');
fprintf(fid,' ========================================================================\n');
fprintf(fid,'%6i %11.2f %9.2f %10.2f %9.2f %10.2f %9.2f\n',[(1:st_elementov)' obtsp_elem]');
fprintf(fid,' ========================================================================\n');


% zapis metode racuna na datoteko IZPIS
fprintf(fid,'\n\n\n                ***************\n');
fprintf(fid,'                 METODA RACUNA\n');  
fprintf(fid,'                ***************\n');
fprintf(fid,' ---------------------------------------------------------\n');
fprintf(fid,'  1 -> materialna nelinearnost in geometrijska linearnost\n\n');
fprintf(fid,'  2 -> materialna in geometrijska nelinearnost\n');
if nacin==1
   fprintf(fid,'         (stopnja Gaussove integracije = %i)\n',stopnja);
end   
if nacin==2
   fprintf(fid,'         (stopnja Lobattove integracije = %i)\n',stopnja);
end
fprintf(fid,'         (stopnja interpolac. polinoma za epsilon= %i)\n',polinom(1));
fprintf(fid,'         (stopnja polinoma za kappo= %i)\n\n',polinom(2));
fprintf(fid,' ---------------------------------------------------------\n');
fprintf(fid,'    IZBRANA METODA: %i\n',racunanje);


disp(sprintf('\n -> Berem matriko temperaturnega polja'))
if razp_T==2
	% prebere matriko 2-D temperaturnega polja jeklenih preènih prerezov
	matrikaHC
end
if razp_T==1
   % prebere matriko 1-D temperaturnega polja
   matrikaHCa
end	


% v primeru, da so strižni prerezi = 0 se raèun izvede brez upoštevanja striga
for i=1:size(G_ploscina,2)
   if G_ploscina(i)==0
      G_ploscina(i)=1e15;
   end
end   

if nacin==1
	% uteži in koordinate Gaussovih integracijskih toèk
	% podatki so zapisani v datoteki Gauss.m
	[wg,xg]=gauss(stopnja);
end
if nacin==2
	% uteži in koordinate Lobattovih integracijskih toèk
	% podatki so zapisani v datoteki Lobatto.m
	[wg,xg]=lobatto(stopnja);
end
% stopnja numeriène integracije za ploskovni integral po preènem prerezu
[wg1,xg1]=gauss(3);


% index -> proste prostostne stopnje
% fix -> fiksne (podprte) prostostne stopnje
index=1:3*(st_vozlisc+st_elementov);
index(3*st_elementov+fix)=[];

% sestavljanje celotne togostne matrike konstrukcije (sestav=1) in raèunanje
%  njene determinante
sestav=0;

for element=1:st_elementov
   
   if sestav==1
      % število interpolacijskih toèk za epsilon oziroma kappo za posamezne elemente
   	if isempty(intersect(element,kap_const))==0
      	% element s konstantnim potekom epsilona in kappe
	      stevec1(element,1)=2;
   	   stevec1(element,2)=1;
	   	stevec1(element,3)=1;
	   else
   	   stevec1(element,1)=sum(polinom)+2;
      	stevec1(element,2)=polinom(1)+1;
	      stevec1(element,3)=polinom(2)+1;
      end
   end %sestav==1 
   
   % stopnje interpolacijskih polinomov za epsilon oziroma kappo za posamezne elemente
   if isempty(intersect(element,kap_const))==0
      % element s konstantnim potekom epsilona in kappe
      stevec(element,1)=0;
      stevec(element,2)=0;
	   stevec(element,3)=0;
   else
      stevec(element,1)=sum(polinom);
      stevec(element,2)=polinom(1);
      stevec(element,3)=polinom(2);
   end
   
   
	% sestavljanje transformacijske matrike elementa
	T0(1,1)=(voz(pq(element,3),2)-voz(pq(element,2),2))/(dolzina(element)/100);
	T0(1,2)=(voz(pq(element,3),3)-voz(pq(element,2),3))/(dolzina(element)/100);
	T0(2,1)=(voz(pq(element,3),3)-voz(pq(element,2),3))/(dolzina(element)/100);
	T0(2,2)=-(voz(pq(element,3),2)-voz(pq(element,2),2))/(dolzina(element)/100);
	T0(3,3)=-1;
   
   T=([T0,zeros(3),zeros(3);zeros(3),T0,zeros(3);zeros(3),zeros(3),T0]);
 
 	% transformacijske matrike elementov - 
   % zapisane v 3D matriki T_elementa  
   T_elementa(:,:,element)=T;
   
   
   if sestav==1
      T1=([eye(stevec1(element,1)+3),zeros(stevec1(element,1)+3,6);...\n
			zeros(3,stevec1(element,1)+3),T0,zeros(3);...\n
         zeros(3,stevec1(element,1)+3),zeros(3),T0]);
   
    	% transformacijske matrike elementov - 
   	% zapisane v 3D matriki T_elementa  
   	T1_elementa(1:stevec1(element,1)+9,1:stevec1(element,1)+9,element)=T1;
   end
   
   % integracijske toèke za ENKRATNI integral
	x_ena(element,:)=dolzina(element)/2*(xg+1);
	% integracijske toèke za DVAKRATNI integral
	x_dva(:,:,element)=dolzina(element)/4*(xg+1)'*(xg+1);
	% integracijske toèke za TRIKRATNI integral
	for i=1:size(xg,2)
   	x_tri(:,:,i,element)=dolzina(element)/8*(xg(i)+1)*((xg+1)'*(xg+1));
	end
   
   % raèun vrednosti Lagrangevih polinomov v integracijskih toèkah integralov (za epsilon)
   for n=1:polinom(1)+1  
	   lagrange_ena(n,:,element)=lagrange(x_ena(element,:),n,polinom(1),dolzina(element));
   	lagrange_dva(:,:,n,element)=lagrange(x_dva(:,:,element),n,polinom(1),dolzina(element));
	end
   
   % raèun vrednosti Lagrangevih polinomov v integracijskih toèkah integralov (za kappo)
   for m=1:polinom(2)+1  
   	lagrange_enak(m,:,element)=lagrange(x_ena(element,:),m,polinom(2),dolzina(element));
	   lagrange_dvak(:,:,m,element)=lagrange(x_dva(:,:,element),m,polinom(2),dolzina(element));
      lagrange_trik(:,:,:,m,element)=...\n
         lagrange(x_tri(:,:,:,element),m,polinom(2),dolzina(element));
	end
end %element=1:st_elementov



% ********************
%  Zaèetne vrednosti:
% ********************
% zaèetne vrednosti posplošenih notranjih in vozlišènih pomikov konstrukcije: 
% na posameznem elementu: R1(0),R2(0),R3(0),u(0),w(0),fi(0),u(L),w(L),fi(L)
% na celotni konstrukciji jih je 3 × (število elementov + število vozlišè)
Ru_konstGKSA=zeros(3*(st_elementov+st_vozlisc),1);
% zaèetne vrednosti posplošenih vozlišènih pomikov posameznega elementa
u_elA=zeros(6,st_elementov);
% zaèetne vrednosti posplošenih notranjih pomikov posameznega elementa: N(0),Q(0),M(0)
R_elA=zeros(3,st_elementov);  
% zaèetne vrednosti epsilonov posameznega elementa: eps(1)...eps(N)
eps_elA=zeros(polinom(1)+1,st_elementov); 
% zaèetne vrednosti parametrov za kappo posameznega elementa: kappa(1)...kappa(M)
kappa_elA=zeros(polinom(2)+1,st_elementov);






% kondenzirane desne strani
f_el=zeros(9,st_elementov);
fde=zeros(9,st_elementov);



% koordinate toèk ploskovnega integrala
z=zeros(size(xg1,2),size(xg1,2),max(lamele),st_prerezov);

% zaèetne vrednosti specifiène spremembe dolžine in psevdoukrivljenosti v integracijskih
%   toèkah enkratnega intregrala vzdolž elementov
eps_stari=zeros(st_elementov,stopnja);
kappa_stari=zeros(st_elementov,stopnja);


% ***********************
%  Leseni preèni prerez
% ***********************
% temperaturna deformacija v integracijskih toèkah lesenega preènega
% prereza
Dt_w=zeros(size(xg1,2),size(xg1,2),max(lamele),max(max(st_prrazd')),st_prerezov);

% mehanske deformacije v integracijskih toèkah lesenega preènega prereza 
D_w1=zeros(size(xg1,2),size(xg1,2),max(lamele),max(max(st_prrazd')),stopnja,st_elementov);
% plastièni del mehanske deformacije v integracijskih toèkah lesenega preènega prereza 
D_wp1=zeros(size(xg1,2),size(xg1,2),max(lamele),max(max(st_prrazd')),stopnja,st_elementov);
% spremenljivka alfa v integracijskih toèkah lesenega preènega prereza
alfaw1=zeros(size(xg1,2),size(xg1,2),max(lamele),max(max(st_prrazd')),stopnja,st_elementov);

% napetosti v integracijskih toèkah lesenega preènega prereza 
sigma_w1=zeros(size(xg1,2),size(xg1,2),max(lamele),max(max(st_prrazd')),stopnja,st_elementov);
% meje plastiènosti v integracijskih toèkah lesenega preènega prereza 
sigmawy1=zeros(size(xg1,2),size(xg1,2),max(lamele),max(max(st_prrazd')),stopnja,st_elementov);
% zaostale napetosti ("back stress") v integracijskih toèkah lesenega preènega prereza 
qu_w1=zeros(size(xg1,2),size(xg1,2),max(lamele),max(max(st_prrazd')),stopnja,st_elementov);

% deformacije krèenja v integracijskih toèkah lesenega preènega prereza 
Ds_w1=zeros(size(xg1,2),size(xg1,2),max(lamele),max(max(st_prrazd')),stopnja,st_elementov);
% deformacije lezenja v integracijskih toèkah lesenega preènega prereza 
Dc_w1=zeros(size(xg1,2),size(xg1,2),max(lamele),max(max(st_prrazd')),stopnja,st_elementov);
% mehano-sorptivne deformacije v integracijskih toèkah lesenega preènega prereza 
Dms_w1=zeros(size(xg1,2),size(xg1,2),max(lamele),max(max(st_prrazd')),stopnja,st_elementov);



% definitnost tangentnih togostnih matrik preènih prerezov elementov (v integracijskih toèkah)
defc_pr=zeros(st_elementov,stopnja);

% determinanta togostne matrike konstrukcije
detK_konst=1;

% najveèje dovoljeno pogojenostno število inverznih matrik, da raèun še poteka
st_pogojenosti=1e10;

   
else % zaèetni inkrement razlièen od 1
	% merjenje raèunskega èasa 
	tm=cputime;
	fid = fopen('izpis.out','at');   
   fseek(fid,0,'eof');
   disp(sprintf('\n  -> Nadaljevanje z %i. èasovnim korakom',(zacetek+1)))
end  %zacetek == 0 



% *****************************
%   INKREMENTNA METODA RAÈUNA
% *****************************
for inkrement=(zacetek+1):(zacetek+st_korakov) % število obtežnih korakov
   
   clear ink_cas1
   Podatki
   if inkrement==1
      ink_cas1=0;
   else   
      if avt==1 % èasovni inkrement se izbira samodejno glede na podane èase
         ink_cas=cas(inkrement)-cas(inkrement-1);
      end   
      ink_cas1=ink_cas;
   end   
   
   stev=0;
   stevb=0; % števec bisekcije


   
   % kazalec delitve velikosti parametra "s"
   divergenca=1;
   
   % kazalec bisekcije
   bisekcija=0;

   
  	% *********************************
   while divergenca==1 | bisekcija==1
  	% *********************************   
    
   if bisekcija==0  

		% zaèetni posplošeni vozlišèni pomiki konstrukcije 
   	Ru_konstGKSA(:,inkrement+1)=Ru_konstGKSA(:,inkrement);
      
      % posplošeni vozlišèni pomiki posameznega elementa
	   u_elA(:,:,inkrement+1)=u_elA(:,:,inkrement);
      
		% posplošeni notranji pomiki posameznega elementa: N(0),Q(0),M(0)
		R_elA(:,:,inkrement+1)=R_elA(:,:,inkrement);  

		% epsiloni posameznega elementa: eps(1)...eps(N)
   	eps_elA(:,:,inkrement+1)=eps_elA(:,:,inkrement);

		% kappe posameznega elementa: kappa(1)...kappa(M)
      kappa_elA(:,:,inkrement+1)=kappa_elA(:,:,inkrement);
      
      % èas pri katerem se izvaja raèun
      if inkrement==1
         casA(inkrement)=0;
      else   
         casA(inkrement)=casA(inkrement-1)+ink_cas1;
      end   
      
      % pomožne vrednosti, ki so bodisi enake vrednostim iz prejšnjega inkrementa 
      %   ali pa vrednostim v novem izhodišènem koraku pri izvajanju bisekcije
      sigma_w1p=sigma_w1(:,:,:,:,:,:,inkrement);
      qu_w1p=qu_w1(:,:,:,:,:,:,inkrement);
      D_w1p=D_w1(:,:,:,:,:,:,inkrement);
      D_wp1p=D_wp1(:,:,:,:,:,:,inkrement);
      Ds_w1p=Ds_w1(:,:,:,:,:,:,inkrement);
      Dc_w1p=Dc_w1(:,:,:,:,:,:,inkrement);
      Dms_w1p=Dms_w1(:,:,:,:,:,:,inkrement);
      Dt_wp=Dt_w(:,:,:,:,:,inkrement);
      
      
   else % bisekcija==1   

		% zaèetni posplošeni vozlišèni pomiki konstrukcije 
      Ru_konstGKSA(:,inkrement+1)=Ru_konstGKSAp;
      
      % posplošeni vozlišèni pomiki posameznega elementa
	   u_elA(:,:,inkrement+1)=u_elAp;

		% posplošeni notranji pomiki posameznega elementa: N(0),Q(0),M(0)
		R_elA(:,:,inkrement+1)=R_elAp;  

		% epsiloni posameznega elementa: eps(1)...eps(N)
   	eps_elA(:,:,inkrement+1)=eps_elAp;

		% kappe posameznega elementa: kappa(1)...kappa(M)
      kappa_elA(:,:,inkrement+1)=kappa_elAp;
      
      % èas pri katerem se izvaja raèun
      casA(inkrement)=casAp+ink_cas1;

   end % bisekcija==1


	for t=1:size(cas,2)-1

		if casA(inkrement) >= cas(t) & casA(inkrement) <= cas(t+1)
         
         % linearna interpolacija temperaturnega polja, razporeditve vlage in obtežnega faktorja
         for pl=1:st_tpolj
            if razp_T==2 % 2-dimenzionalno temperaturno polje prereza 
              	for pod=1:size(T_HC,2)
     	         	T_HC1{inkrement,pod,pl}=T_HC{t,pod,pl}+...\n
    		         	(T_HC{t+1,pod,pl}-T_HC{t,pod,pl})*((casA(inkrement)-cas(t))/(cas(t+1)-cas(t)));
                    W_HC1{inkrement,pod,pl}=W_HC{t,pod,pl}+...\n
    		         	(W_HC{t+1,pod,pl}-W_HC{t,pod,pl})*((casA(inkrement)-cas(t))/(cas(t+1)-cas(t)));
	            end
            else % konstantna temperatura prereza (razp_T==1)
               T_HCa1(inkrement,pl)=T_HCa(t,pl)+...\n
       	         (T_HCa(t+1,pl)-T_HCa(t,pl))*((casA(inkrement)-cas(t))/(cas(t+1)-cas(t)));
            end   
	      end %pl=1:size(T_HC,3)
   
   	   obt_faktor1(inkrement+1)=obt_faktor(t)+(obt_faktor(t+1)-...\n
      	   obt_faktor(t))*((casA(inkrement)-cas(t))/(cas(t+1)-cas(t)));
         
         break % prekine se zanka t=1:size(cas,2)-1
      end   
   end %t=1:size(cas,2)-1

   %inkrement posplošenih vozlišènih pomikov: 
	% na posameznem elementu: du(0),dw(0),dfi(0),du(L),dw(L),dfi(L)
	% na celotni konstrukciji jih je 3 × število vozlišè
   du_konstGKS=zeros(3*st_vozlisc,1);
   

   % indikator, kdaj zaènemo sestavljati celotno togostno matriko konstrukcije (raèun detK)
   zacni=0;
  
   
   disp(sprintf('\n\n\n  --------------------------------------------------'))
   disp(sprintf('   ->  Raèunski korak = %i  (t= %.2f h)',inkrement,casA(inkrement)/3600))
   disp('  --------------------------------------------------')
   
   
   % *****************************************************************************************
	%  RAÈUN TEMPERATURNEGA POLJA TER KARAKTERISTIK MATERIALOV V INTEG. TOÈKAH PREÈNEGA PREREZA
	% *****************************************************************************************
   for p=1:st_prerezov
      % temperaturno polje obravnavanega preènega prereza
		if isempty(tpr)
   		polje=1;
		else   
			polje=tpr(p,2);
      end
      
      if razp_T==1
         T_HC1(inkrement,:,polje)=0;
         W_HC1(inkrement,:,polje)=0;
         dy(polje,:)=0;
         dz(polje,:)=0;
      end
      if razp_T==2
         T_HCa1(inkrement,polje)=0;
      end   
      
      % raèun se izvede le enkrat za vsak posamezen raèunski korak
      % zaradi simetriène razporeditve temperature se obravnava le polovica preènega prereza
      [z(:,:,1:lamele(p),p),Dt_w(:,:,1:lamele(p),1:max(st_prrazd(p,:)),p,inkrement+1),...\n
      Et_T(:,:,1:lamele(p),:,p),Ec_T(:,:,1:lamele(p),:,p),...\n
      ft_T(:,:,1:lamele(p),:,p),fc_T(:,:,1:lamele(p),:,p),...\n
      temp_w(:,:,:,:,p,inkrement+1),moist_w(:,:,:,:,p,inkrement+1)]=...\n
   		temperaturaHC(st_lamel,p,st_prerezov,st_podprerezov(p),visina(p,:),prrazd(p,:),zT_lamele_zg,...\n
        zT_lamele_sp,b_lamele,zT(p),wg1,xg1,dy(polje,:),dz(polje,:),T_HC1(inkrement,:,polje),W_HC1(inkrement,:,polje),...\n
        T_HCa1(inkrement,polje),E_t0(p),E_c0(p),f_t0(p),f_c0,razp_T);
    
   end %p=1:st_prerezov
   
   % sprememba temperaturne deformacije v posameznem èasovnem inkrementu
   dDt_w=Dt_w(:,:,:,:,:,inkrement+1)-Dt_wp;
   % sprememba vlage v posameznem èasovnem inkrementu
   dmoist_w=moist_w(:,:,:,:,:,inkrement+1)-moist_w(:,:,:,:,:,inkrement);

   
   if inkrement==1
      EtT=Et_T;
      EcT=Ec_T;
      dmoist_w=zeros(size(xg1,2),size(xg1,2),max(lamele),max(max(st_prrazd')),st_prerezov);
   end
      
   
   disp(sprintf('\n ------------------------------------------------------------------------'));
   if sestav==1
      disp(sprintf('%8s %15s %15s %14s %17s','korak','NORMA delta_u','NORMA f_konst','maxNORMA_el','det(K1_konst)'));
	else   
      disp(sprintf('%8s %15s %15s %14s %16s','korak','NORMA delta_u','NORMA f_konst','maxNORMA_el','det(K_konst)'));
   end   
   disp(' ------------------------------------------------------------------------')
   

   
   
   % ***************************************************************
   %  Raèun inkrementa spremenljivk z Newtonovo iteracijsko metodo 
   % ***************************************************************
   for iter=1:25 %20
      
      % kazalec delitve velikosti parametra "s"
      divergenca=0;
      
      
      % zaèetna togostna matrika je enaka nièelni matriki
      K_konst=zeros(3*(st_elementov+st_vozlisc));
      % uporaba razpršenih matrik
      %K_konst=sparse(3*(st_elementov+st_vozlisc),3*(st_elementov+st_vozlisc));
      
      if sestav==1
	      % celotna togostna matrika konstrukcije
         K1_konst=zeros(sum(stevec1(1:st_elementov,1))+3*(st_elementov+st_vozlisc));
         % uporaba razpršenih matrik
         %K1_konst=sparse(sum(stevec1(1:st_elementov,1))+3*(st_elementov+st_vozlisc),...\n
         %   sum(stevec1(1:st_elementov,1))+3*(st_elementov+st_vozlisc));
      end
      
	
		% zaèetni vektor desnih strani je enak nièelnemu vektorju
      f_konst=zeros(3*(st_elementov+st_vozlisc),1);
      

      % element s konstantnim potekom epsilona in kappe
      if isempty(intersect(element,kap_const))==0
	      K_1o=zeros(2,11-1,st_elementov);
      else
      	K_1o=zeros(sum(polinom)+2,sum(polinom)+11-1,st_elementov);
      end
      
      % vrednosti mej plastiènega teèenja v integracijskih toèkah prereza
      %  izraèunajo se samo v prvem koraku posameznega inkrementa
      if iter==1
         sigwy1=sigmawy1;
      else
         sigwy1=sigma_wy1;
      end   
	   
      for element=1:st_elementov

         % deformacijske,konstitucijske in geometrijske kolièine obravnavanega elementa
         prirast(:,element)=[eps_elA(:,element,inkrement+1);kappa_elA(:,element,inkrement+1);...\n
         	R_elA(:,element,inkrement+1);u_elA(:,element,inkrement+1)];
         
	      % preèni prerez obravnavanega elementa
	      if isempty(pr)
	         prerez=1;
	      else   
				prerez=pr(element,2);
	      end
	      
			% *************************************************************
			%  Raèun togostne matrike in obtežnega vektorja elementa v LKS
         % *************************************************************
         % racunanje=1 -> geometrijska linearnost in materialna nelinearnost
         % racunanje=2 -> geometrijska in materialna nelinearnost 
         
         if racunanje==1
            pause
         end   
        
        if racunanje==2
            
            % element s konstantnim potekom epsilona in kappe
            if isempty(intersect(element,kap_const))==0
             
               
            else % element, pri katerem sta epsilon in kappa interp. z Lagrangevimi polinomi
            
                % funkcija K_GMN je zapisana v datoteki K_GMN.m 
	            [K_ee(1:stevec(element,2)+1,1:stevec(element,2)+1,element),...\n
   	            K_ek(1:stevec(element,2)+1,1:stevec(element,3)+1,element),...\n
      	        K_eR(1:stevec(element,2)+1,:,element),K_eu(1:stevec(element,2)+1,:,element),...\n
         	    K_ke(1:stevec(element,3)+1,1:stevec(element,2)+1,element),...\n
            	K_kk(1:stevec(element,3)+1,1:stevec(element,3)+1,element),...\n
	            K_kR(1:stevec(element,3)+1,:,element),K_ku(1:stevec(element,3)+1,:,element),...\n
   	            K_Re(:,1:stevec(element,2)+1,element),K_Rk(:,1:stevec(element,3)+1,element),...\n
      	        K_RR(:,:,element),K_Ru(:,:,element),K_ue(:,1:stevec(element,2)+1,element),...\n
         	    K_uk(:,1:stevec(element,3)+1,element),K_uR(:,:,element),K_uu(:,:,element),...\n
            	g_e(1:stevec(element,2)+1,element),g_k(1:stevec(element,3)+1,element),...\n
                f_R(:,element),f_u(:,element),alfa_w1(:,:,:,:,:,element),...\n
                D_w1(:,:,:,:,:,element,inkrement+1),D_wp1(:,:,:,:,:,element,inkrement+1),...\n
   	            sigma_w1(:,:,:,:,:,element,inkrement+1),sigma_wy1(:,:,:,:,:,element),...\n
                qu_w1(:,:,:,:,:,element,inkrement+1),Ds_w1(:,:,:,:,:,element,inkrement+1),...\n
                Dc_w1(:,:,:,:,:,element,inkrement+1),Dms_w1(:,:,:,:,:,element,inkrement+1),...\n
                divergenca,eps_ena(element,:),kappa_ena(element,:),Nc1(element,:,inkrement+1),...\n
                Mc1(element,:,inkrement+1),defc_pr(element,:,inkrement+1),...\n
                c11_ena(element,:,inkrement),c12_ena(element,:,inkrement),c22_ena(element,:,inkrement)]=...\n        
      			    K_GMNles(element,st_prerezov,prerez,dolzina(element),polinom,stopnja,...\n
                x_ena(element,:),x_dva(:,:,element),x_tri(:,:,:,element),...\n
                lagrange_ena(:,:,element),lagrange_dva(:,:,:,element),...\n
                lagrange_enak(:,:,element),lagrange_dvak(:,:,:,element),...\n
                lagrange_trik(:,:,:,:,element),xg,wg,xg1,wg1,eps_stari(element,:),...\n
                kappa_stari(element,:),prirast(:,element)',st_lamel,prrazd(prerez,:),b_lamele,...\n
                zT_lamele_sp,zT_lamele_zg,z(:,:,:,prerez),alfaw1(:,:,:,:,:,element),...\n
                D_w1p(:,:,:,:,:,element),D_wp1p(:,:,:,:,:,element),sigma_w1p(:,:,:,:,:,element),...\n
                sigwy1(:,:,:,:,:,element),qu_w1p(:,:,:,:,:,element),EtT(:,:,:,:,prerez),...\n
                Et_T(:,:,:,:,prerez),EcT(:,:,:,:,prerez),Ec_T(:,:,:,:,prerez),temp_w(:,:,:,:,prerez,inkrement+1),...\n
                dmoist_w(:,:,:,:,prerez),dDt_w(:,:,:,:,prerez),krcenje_w,lezenje_w,...\n
                mod_cw,mehso_w,Ds_w1p(:,:,:,:,:,element),Dc_w1p(:,:,:,:,:,element),...\n
                Dms_w1p(:,:,:,:,:,element),ft_T(:,:,:,:,prerez),fc_T(:,:,:,:,prerez),...\n
                G_modul(prerez),G_ploscina(prerez),obt_faktor1(inkrement+1),...\n
                obtn_elem(element,:),obtsp_elem(element,:),diagram(prerez),inkrement,iter,fid,casA);
            
           		if divergenca==1
            		break % prekine se zanka element=1:st_elementov
		         end
   	      end % element, pri katerem sta epsilon in kappa interp. z lagrangevim polinomom
            
         end %racunanje==2
         
         
         
         % ***********************************************
			%  Togostna matrika konstrukcije in desne strani
			% ***********************************************
         if sestav==1 & zacni==1
           
            % transformacija togostne matrike elementa iz LKS v GKS
            K1_uGKS=T1_elementa(1:stevec1(element,1)+9,1:stevec1(element,1)+9,element)*...
            	([K_ee(1:stevec(element,2)+1,1:stevec(element,2)+1,element),...\n
					K_ek(1:stevec(element,2)+1,1:stevec(element,3)+1,element),...\n
      	      K_eR(1:stevec(element,2)+1,:,element),K_eu(1:stevec(element,2)+1,:,element);...\n
   	         K_ke(1:stevec(element,3)+1,1:stevec(element,2)+1,element),...\n
	            K_kk(1:stevec(element,3)+1,1:stevec(element,3)+1,element),...\n
            	K_kR(1:stevec(element,3)+1,:,element),K_ku(1:stevec(element,3)+1,:,element);...\n
         	   K_Re(:,1:stevec(element,2)+1,element),K_Rk(:,1:stevec(element,3)+1,element),...\n
      	      K_RR(:,:,element) K_Ru(:,:,element);...\n
   	         K_ue(:,1:stevec(element,2)+1,element),K_uk(:,1:stevec(element,3)+1,element),...\n
	            K_uR(:,:,element) K_uu(:,:,element)])*...\n
         	   T1_elementa(1:stevec1(element,1)+9,1:stevec1(element,1)+9,element)';
         
           
				% ********************************************
				%  Sestavljanje togostne matrike konstrukcije
   	      % ********************************************         
      	   j=1:3;
         	r=1:stevec1(element,2);
	         s=1:stevec1(element,3);
         
   	      % K_ee
      	   K1_konst(sum(stevec1(1:element-1,2))+r,sum(stevec1(1:element-1,2))+r)=...\n
               K1_konst(sum(stevec1(1:element-1,2))+r,sum(stevec1(1:element-1,2))+r)+...\n
               K1_uGKS(r,r);
	         % K_ek
   	      K1_konst(sum(stevec1(1:element-1,2))+r,sum(stevec1(1:st_elementov,2))+...\n
      	      sum(stevec1(1:element-1,3))+s)=...\n
            	K1_konst(sum(stevec1(1:element-1,2))+r,sum(stevec1(1:st_elementov,2))+...\n
            	sum(stevec1(1:element-1,3))+s)+K1_uGKS(r,max(r)+s);
         	% K_eR
	         K1_konst(sum(stevec1(1:element-1,2))+r,sum(stevec1(1:st_elementov,1))+...\n
   	         3*(element-1)+j)=...\n
      	      K1_konst(sum(stevec1(1:element-1,2))+r,sum(stevec1(1:st_elementov,1))+...\n
         	   3*(element-1)+j)+K1_uGKS(r,max(r)+max(s)+j);
         	% K_eui
	         K1_konst(sum(stevec1(1:element-1,2))+r,sum(stevec1(1:st_elementov,1))+...\n
   	         3*(st_elementov+pq(element,2)-1)+j)=...\n
      	      K1_konst(sum(stevec1(1:element-1,2))+r,sum(stevec1(1:st_elementov,1))+...\n
         	   3*(st_elementov+pq(element,2)-1)+j)+K1_uGKS(r,max(r)+max(s)+j+3);
	         % K_euj
   	      K1_konst(sum(stevec1(1:element-1,2))+r,sum(stevec1(1:st_elementov,1))+...\n
      	      3*(st_elementov+pq(element,3)-1)+j)=...\n
         	   K1_konst(sum(stevec1(1:element-1,2))+r,sum(stevec1(1:st_elementov,1))+...\n
            	3*(st_elementov+pq(element,3)-1)+j)+K1_uGKS(r,max(r)+max(s)+j+6);
         
         	% K_ke
	         K1_konst(sum(stevec1(1:st_elementov,2))+sum(stevec1(1:element-1,3))+s,...\n
   	         sum(stevec1(1:element-1,2))+r)=...\n
      	      K1_konst(sum(stevec1(1:st_elementov,2))+sum(stevec1(1:element-1,3))+s,...\n
         	   sum(stevec1(1:element-1,2))+r)+K1_uGKS(max(r)+s,r);
         	% K_kk
	         K1_konst(sum(stevec1(1:st_elementov,2))+sum(stevec1(1:element-1,3))+s,...\n
   	         sum(stevec1(1:st_elementov,2))+sum(stevec1(1:element-1,3))+s)=...\n
      	      K1_konst(sum(stevec1(1:st_elementov,2))+sum(stevec1(1:element-1,3))+s,...\n
               sum(stevec1(1:st_elementov,2))+sum(stevec1(1:element-1,3))+s)+...\n
               K1_uGKS(max(r)+s,max(r)+s);
	         % K_kR
   	      K1_konst(sum(stevec1(1:st_elementov,2))+sum(stevec1(1:element-1,3))+s,...\n
      	      sum(stevec1(1:st_elementov,1))+3*(element-1)+j)=...\n
         	   K1_konst(sum(stevec1(1:st_elementov,2))+sum(stevec1(1:element-1,3))+s,...\n
               sum(stevec1(1:st_elementov,1))+3*(element-1)+j)+...\n
               K1_uGKS(max(r)+s,max(r)+max(s)+j);
         	% K_kui
	         K1_konst(sum(stevec1(1:st_elementov,2))+sum(stevec1(1:element-1,3))+s,...\n
   	         sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,2)-1)+j)=...\n
      	      K1_konst(sum(stevec1(1:st_elementov,2))+sum(stevec1(1:element-1,3))+s,...\n
         	   sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,2)-1)+j)+...\n
            	K1_uGKS(max(r)+s,max(r)+max(s)+j+3);
         	% K_kuj
	         K1_konst(sum(stevec1(1:st_elementov,2))+sum(stevec1(1:element-1,3))+s,...\n
   	         sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,3)-1)+j)=...\n
      	      K1_konst(sum(stevec1(1:st_elementov,2))+sum(stevec1(1:element-1,3))+s,...\n
         	   sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,3)-1)+j)+...\n
            	K1_uGKS(max(r)+s,max(r)+max(s)+j+6);
         
         	% K_Re
	         K1_konst(sum(stevec1(1:st_elementov,1))+3*(element-1)+j,...\n
   	         sum(stevec1(1:element-1,2))+r)=...\n
      	      K1_konst(sum(stevec1(1:st_elementov,1))+3*(element-1)+j,...\n
         	   sum(stevec1(1:element-1,2))+r)+K1_uGKS(max(r)+max(s)+j,r);
         	% K_Rk
	         K1_konst(sum(stevec1(1:st_elementov,1))+3*(element-1)+j,...\n
   	         sum(stevec1(1:st_elementov,2))+sum(stevec1(1:element-1,3))+s)=...\n
      	      K1_konst(sum(stevec1(1:st_elementov,1))+3*(element-1)+j,...\n
               sum(stevec1(1:st_elementov,2))+sum(stevec1(1:element-1,3))+s)+...\n
               K1_uGKS(max(r)+max(s)+j,max(r)+s);
         	% K_RR
	         K1_konst(sum(stevec1(1:st_elementov,1))+3*(element-1)+j,...\n
   	         sum(stevec1(1:st_elementov,1))+3*(element-1)+j)=...\n
      	      K1_konst(sum(stevec1(1:st_elementov,1))+3*(element-1)+j,...\n
               sum(stevec1(1:st_elementov,1))+3*(element-1)+j)+...\n
               K1_uGKS(max(r)+max(s)+j,max(r)+max(s)+j);
	         % K_Rui
   	      K1_konst(sum(stevec1(1:st_elementov,1))+3*(element-1)+j,...\n
      	      sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,2)-1)+j)=...\n
         	   K1_konst(sum(stevec1(1:st_elementov,1))+3*(element-1)+j,...\n
            	sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,2)-1)+j)+...\n
            	K1_uGKS(max(r)+max(s)+j,max(r)+max(s)+j+3);
	         % K_Ruj
   	      K1_konst(sum(stevec1(1:st_elementov,1))+3*(element-1)+j,...\n
      	      sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,3)-1)+j)=...\n
         	   K1_konst(sum(stevec1(1:st_elementov,1))+3*(element-1)+j,...\n
            	sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,3)-1)+j)+...\n
	            K1_uGKS(max(r)+max(s)+j,max(r)+max(s)+j+6);
         
         	% K_uie
	         K1_konst(sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,2)-1)+j,...\n
   	         sum(stevec1(1:element-1,2))+r)=...\n
      	      K1_konst(sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,2)-1)+j,...\n
         	   sum(stevec1(1:element-1,2))+r)+K1_uGKS(max(r)+max(s)+j+3,r);
         	% K_uik
	         K1_konst(sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,2)-1)+j,...\n
   	         sum(stevec1(1:st_elementov,2))+sum(stevec1(1:element-1,3))+s)=...\n
      	      K1_konst(sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,2)-1)+j,...\n
               sum(stevec1(1:st_elementov,2))+sum(stevec1(1:element-1,3))+s)+...\n
               K1_uGKS(max(r)+max(s)+j+3,max(r)+s);
	         % K_uiR
   	      K1_konst(sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,2)-1)+j,...\n
      	      sum(stevec1(1:st_elementov,1))+3*(element-1)+j)=...\n
         	   K1_konst(sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,2)-1)+j,...\n
               sum(stevec1(1:st_elementov,1))+3*(element-1)+j)+...\n
               K1_uGKS(max(r)+max(s)+j+3,max(r)+max(s)+j);
         	% K_uiui
	         K1_konst(sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,2)-1)+j,...\n
   	         sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,2)-1)+j)=...\n
      	      K1_konst(sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,2)-1)+j,...\n
         	   sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,2)-1)+j)+...\n
            	K1_uGKS(max(r)+max(s)+j+3,max(r)+max(s)+j+3);
         	% K_uiuj
	         K1_konst(sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,2)-1)+j,...\n
   	         sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,3)-1)+j)=...\n
      	      K1_konst(sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,2)-1)+j,...\n
         	   sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,3)-1)+j)+...\n
            	K1_uGKS(max(r)+max(s)+j+3,max(r)+max(s)+j+6);
         
         	% K_uje
	         K1_konst(sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,3)-1)+j,...\n
   	         sum(stevec1(1:element-1,2))+r)=...\n
      	      K1_konst(sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,3)-1)+j,...\n
         	   sum(stevec1(1:element-1,2))+r)+K1_uGKS(max(r)+max(s)+j+6,r);
	         % K_ujk
   	      K1_konst(sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,3)-1)+j,...\n
      	      sum(stevec1(1:st_elementov,2))+sum(stevec1(1:element-1,3))+s)=...\n
         	   K1_konst(sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,3)-1)+j,...\n
               sum(stevec1(1:st_elementov,2))+sum(stevec1(1:element-1,3))+s)+...\n
               K1_uGKS(max(r)+max(s)+j+6,max(r)+s);
         	% K_ujR
	         K1_konst(sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,3)-1)+j,...\n
   	         sum(stevec1(1:st_elementov,1))+3*(element-1)+j)=...\n
      	      K1_konst(sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,3)-1)+j,...\n
               sum(stevec1(1:st_elementov,1))+3*(element-1)+j)+...\n
               K1_uGKS(max(r)+max(s)+j+6,max(r)+max(s)+j);
         	% K_ujui
	         K1_konst(sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,3)-1)+j,...\n
   	         sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,2)-1)+j)=...\n
      	      K1_konst(sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,3)-1)+j,...\n
         	   sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,2)-1)+j)+...\n
            	K1_uGKS(max(r)+max(s)+j+6,max(r)+max(s)+j+3);
         	% K_ujuj
	         K1_konst(sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,3)-1)+j,...\n
   	         sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,3)-1)+j)=...\n
      	      K1_konst(sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,3)-1)+j,...\n
         	   sum(stevec1(1:st_elementov,1))+3*(st_elementov+pq(element,3)-1)+j)+...\n
            	K1_uGKS(max(r)+max(s)+j+6,max(r)+max(s)+j+6);

 			end % sestav==1 & zacni==1
         
         % **************************************************
			%  Kondenzacija tangentne togostne matrike elementa
			% **************************************************
         % kondenzacija vsakega epsilona in kappe posebej
      
        	K_oo=...\n
         	([K_ee(1:stevec(element,2)+1,1:stevec(element,2)+1,element),...\n
				K_ek(1:stevec(element,2)+1,1:stevec(element,3)+1,element),...\n
            K_eR(1:stevec(element,2)+1,:,element),K_eu(1:stevec(element,2)+1,:,element);...\n
            K_ke(1:stevec(element,3)+1,1:stevec(element,2)+1,element),...\n
            K_kk(1:stevec(element,3)+1,1:stevec(element,3)+1,element),...\n
            K_kR(1:stevec(element,3)+1,:,element),K_ku(1:stevec(element,3)+1,:,element);...\n
            K_Re(:,1:stevec(element,2)+1,element),K_Rk(:,1:stevec(element,3)+1,element),...\n
            K_RR(:,:,element) K_Ru(:,:,element);...\n
            K_ue(:,1:stevec(element,2)+1,element),K_uk(:,1:stevec(element,3)+1,element),...\n
            K_uR(:,:,element) K_uu(:,:,element)]);
        

			f_o=([-g_e(1:stevec(element,2)+1,element);-g_k(1:stevec(element,3)+1,element);...\n
         	-f_R(:,element);-f_u(:,element)]);
      
      
      
      	for m=1:stevec(element,1)+2
             
            K_11(element,m)=K_oo(1,1);
            K_1o(m,1:stevec(element,1)+11-m,element)=K_oo(1,2:size(K_oo,2));
            f1(element,m)=f_o(1);
               
            if K_oo(1,1) == 0  &  inkrement > 1
	     			divergenca=1; 
	  	   		break  % prekine se zanka  m=1:stevec(element,1)+5
            end
            
            KK_oo=K_oo(2:size(K_oo,2),2:size(K_oo,2))-...\n
               K_oo(2:size(K_oo,2),1)*(1/K_oo(1,1))*K_oo(1,2:size(K_oo,2));
            ff_o=f_o(2:size(f_o,1),1)-...\n
               K_oo(2:size(K_oo,2),1)*(1/K_oo(1,1))*f_o(1);
              
            clear K_oo do_dlambda f_o
            K_oo=KK_oo;
            f_o=ff_o;
               
            if m==stevec(element,1)+2
               % togostna matrika elementa je dimenzije [9×9]
               K_el(:,:,element)=K_oo;
               f_el(:,element)=f_o;
               clear K_oo f_o
            end   
            clear KK_oo ff_o 
            
         end % m=1:stevec(element,1)+2 
         
         if divergenca == 1
            break  % prekine se element=1:st_elementov
         end
         

         
			% kondenzacija togostne matrike elementa zaradi sprostitve izbranih prost. stopenj
			% funkcija je zapisana v datoteki "kondenzac.m"
         if sprost(element,:)==0 
            K_el1(:,:,element)=K_el(:,:,element);
            if isempty(spr)==0
               f_b(:,element)=zeros(max(s_p),1);
            end   
         else 
            pause
         	K_el1(:,:,element)=K_el(:,:,element);
            fde(:,element)=f_el(:,element);
            [K_el(:,:,element),f_el(:,element),f_b(:,element)]=...\n
               kondenzac(K_el1(:,:,element),fde(:,element),sprost(element,:));
         end


         
      	% *********************************
			%  Togostna matrika elementa v GKS
      	% *********************************
      	% transformacija togostna matrike elementa 
         K_uGKS=T_elementa(:,:,element)*K_el(:,:,element)*T_elementa(:,:,element)';
         
	
			% *******************************************
			%  Sestavljanje togostne matrike konstrukcije
			% *******************************************
			j=1:3;
         k=1:3;

			% združevanje po posplošenih notranjih pomikih
			K_konst(3*(element-1)+j,3*(element-1)+k)=...\n
            K_konst(3*(element-1)+j,3*(element-1)+k)+K_uGKS(j,k);
         K_konst(3*(element-1)+j,3*(st_elementov+pq(element,2)-1)+k)=...\n
            K_konst(3*(element-1)+j,3*(st_elementov+pq(element,2)-1)+k)+K_uGKS(j,k+3);
         K_konst(3*(element-1)+j,3*(st_elementov+pq(element,3)-1)+k)=...\n
            K_konst(3*(element-1)+j,3*(st_elementov+pq(element,3)-1)+k)+K_uGKS(j,k+6);
         K_konst(3*(st_elementov+pq(element,2)-1)+j,3*(element-1)+k)=...\n
            K_konst(3*(st_elementov+pq(element,2)-1)+j,3*(element-1)+k)+K_uGKS(j+3,k);
         K_konst(3*(st_elementov+pq(element,3)-1)+j,3*(element-1)+k)=...\n
            K_konst(3*(st_elementov+pq(element,3)-1)+j,3*(element-1)+k)+K_uGKS(j+6,k);

			% združevanje po posplošenih vozlišènih pomikih
			K_konst(3*(st_elementov+pq(element,2)-1)+j,3*(st_elementov+pq(element,2)-1)+k)=...\n
            K_konst(3*(st_elementov+pq(element,2)-1)+j,3*(st_elementov+pq(element,2)-1)+k)+...\n
            K_uGKS(j+3,k+3);
		   K_konst(3*(st_elementov+pq(element,2)-1)+j,3*(st_elementov+pq(element,3)-1)+k)=...\n
            K_konst(3*(st_elementov+pq(element,2)-1)+j,3*(st_elementov+pq(element,3)-1)+k)+...\n
            K_uGKS(j+3,k+6);
		   K_konst(3*(st_elementov+pq(element,3)-1)+j,3*(st_elementov+pq(element,2)-1)+k)=...\n
            K_konst(3*(st_elementov+pq(element,3)-1)+j,3*(st_elementov+pq(element,2)-1)+k)+...\n
            K_uGKS(j+6,k+3);
		   K_konst(3*(st_elementov+pq(element,3)-1)+j,3*(st_elementov+pq(element,3)-1)+k)=...\n
            K_konst(3*(st_elementov+pq(element,3)-1)+j,3*(st_elementov+pq(element,3)-1)+k)+...\n
            K_uGKS(j+6,k+6);
         
            
			% *******************************
			%  Obtežni vektor elementa v GKS
      	% *******************************
      	f_elGKS=T_elementa(:,:,element)*f_el(:,element);
	
			% **********************************************
			%  Sestavljanje obtežnega vektorja konstrukcije
			% **********************************************
         j=1:3;
         
         % združevanje po posplošenih notranjih pomikih
      	f_konst(3*(element-1)+j)=...\n
            f_konst(3*(element-1)+j)+f_elGKS(j);
         
      	% združevanje po posplošenih vozlišènih pomikih
         f_konst(3*(st_elementov+pq(element,2)-1)+j)=...\n
            f_konst(3*(st_elementov+pq(element,2)-1)+j)+f_elGKS(j+3);
         f_konst(3*(st_elementov+pq(element,3)-1)+j)=...\n
            f_konst(3*(st_elementov+pq(element,3)-1)+j)+f_elGKS(j+6);
         
      end %element=1:st_elementov
      
      %if divergenca==1
      %   break %prekine se zanka iter=1:20
      %end   
      
      
      % nespremenljiva in spremenljiva vozlišèna obtežba konstrukcije
      obtn_voz1=obtn_voz';
      obtsp_voz1=obtsp_voz';
      % spreminjanje enot pri vozlišèni momentni obtežbi
      obtn_voz1(3,:)=obtn_voz1(3,:)*100;
      obtsp_voz1(3,:)=obtsp_voz1(3,:)*100;
      
		% skupni obtežni vektor
      f_konst(3*(st_elementov)+1:3*(st_elementov+st_vozlisc))=...\n
         f_konst(3*(st_elementov)+1:3*(st_elementov+st_vozlisc))+...\n
         obtn_voz1(:)+obt_faktor1(inkrement+1)*obtsp_voz1(:);
      
      

		% ********************************************************************
		%  Prirastek posplošenih vozlišènih in notranjih pomikov konstrukcije 
		% ********************************************************************
	   
	   % èrtamo prostostne stopnje, kjer so vozlišèni pomiki predpisani
		% fix -> fiksne(podprte) prostostne stopnje
		K_konst(:,3*st_elementov+fix)=[];
		K_konst(3*st_elementov+fix,:)=[];
      f_konst(3*st_elementov+fix)=[];
      
      if sestav==1 & zacni==1
         K1_konst(:,sum(stevec1(1:st_elementov,1))+3*st_elementov+fix)=[];
			K1_konst(sum(stevec1(1:st_elementov,1))+3*st_elementov+fix,:)=[];
      end
      
      % determinanta togostne matrike konstrukcije
      if sestav==1 & zacni==1
         detK_konst(inkrement+1)=det(K1_konst);
         %detK_konst(inkrement+1)=2.66-casA(inkrement);
      else
         detK_konst(inkrement+1)=det(K_konst);
      	if abs(detK_konst(inkrement+1)) < 1e-50 & iter > 1 & inkrement > 1
         	divergenca=1;
	      end 
      end
         
      
      deltaRu_konst=zeros(3*(st_elementov+st_vozlisc),1);
      % prirastki posplošenih notranjih in vozlišènih pomikov konstrukcije
      deltaRu_konst(index)=K_konst\f_konst;
          
      
	   % Raèun prirastkov posplošenih notranjih pomikov dN(0),dQ(0),dM(0) posameznega elementa 
		% Sestavljanje vektorja {dNi(0),dQi(0),dMi(0)} na nivoju konstrukcije
	   for element=1:st_elementov
	      
			% prirastki posplošenih notranjih in vozlišènih pomikov posameznega elementa v GKS
         deltaRu_el(:,element)=...\n
            deltaRu_konst([(3*(element-1)+1):(3*element),...\n
            ((st_elementov+pq(element,2)-1)*3+1):((st_elementov+pq(element,2))*3),...\n
            ((st_elementov+pq(element,3)-1)*3+1):((st_elementov+pq(element,3))*3)]);
            
         % transformacija v LKS
         deltaRu_el(:,element)=T_elementa(:,:,element)'*deltaRu_el(:,element);
            
         % prirastki posplošenih notranjih pomikov
         deltaR_el(:,element)=deltaRu_el(1:3,element);
            
         % prirastki posplošenih vozlišènih pomikov
         deltau_el(:,element)=deltaRu_el(4:9,element);
         
         
         % raèun prirastkov posplošenih vozlišènih pomikov 
         %   na mestu sprostitve prostostne stopnje
	   	if sum(sprost(element,:)) > 0
   	      delu_el=deltau_el(:,element);
	   	  	delu_el(sprost(element,:))=[];
            in=1:size(K_el1(:,:,element),1);
      		in(:,sprost(element,:))=[];
            
            %deltau_el=inv(K_bb)*f_b-inv(K_bb)*K_ba*delu_el-...\n
           	%               inv(K_bb)*flambda_b*delta_lambda;
			   deltau_el(sprost(element,:),element)=...\n
   	     		inv(K_el1(sprost(element,:),sprost(element,:),element))*f_b(:,element)-...\n
      	     	inv(K_el1(sprost(element,:),sprost(element,:),element))*...\n
               (K_el1(sprost(element,:),in,element))*delu_el-...\n
           	   inv(K_el1(sprost(element,:),sprost(element,:),element))*...\n
               flambda_b(:,element)*delta_lambda;
	      end %sum(sprost(element,:)) > 0
         
         
         % kondenzacija   
         delta_all=deltaRu_el(:,element);
         for m=1:stevec(element,1)+2
               delta1=1/K_11(element,stevec(element,1)+2+1-m)*...\n
                  (f1(element,stevec(element,1)+2+1-m)-...\n
                  K_1o(stevec(element,1)+2+1-m,...\n
                  1:stevec(element,1)+11-(stevec(element,1)+2+1-m),element)*delta_all);
               delta_all=([delta1;delta_all]);

                  
					if m==(stevec(element,1)+2)
                  delta1e_el=delta_all(1:stevec(element,2)+1);
                  delta1k_el=delta_all(stevec(element,2)+2:stevec(element,1)+2);
                     
                  if isempty(intersect(element,kap_const))==0
                     % element s konstantnim potekom epsilona in kappe
                     deltae_el(1:polinom(1)+1,element)=...\n
                        ones(polinom(1)+1,1).*delta1e_el;
                     deltak_el(1:polinom(2)+1,element)=...\n
                        ones(polinom(2)+1,1).*delta1k_el;
                  else
                     deltae_el(:,element)=delta1e_el;
                     deltak_el(:,element)=delta1k_el;
                  end   
                  clear delta_all delta1e_el  delta1k_el
               end   
               
            end % m=1:stevec(element,1)+2
         
      end %element=1:st_elementov
      
     
      % popravljeni posplošeni notranji in vozlišèni pomiki konstrukcije 
      Ru_konstGKSA(:,inkrement+1)=Ru_konstGKSA(:,inkrement+1)+deltaRu_konst;
      
      % popravljene vrednosti inkrementa posplošenih vozlišènih pomikov konstrukcije
      du_konstGKS=du_konstGKS+deltaRu_konst(3*st_elementov+1:3*(st_elementov+st_vozlisc));
      
      % popravljeni posplošeni vozlišèni pomiki posameznega elementa
      u_elA(:,:,inkrement+1)=u_elA(:,:,inkrement+1)+deltau_el;

		% popravljeni posplošeni notranji pomiki posameznega elementa: N(0),Q(0),M(0)
		R_elA(:,:,inkrement+1)=R_elA(:,:,inkrement+1)+deltaR_el;  

		% popravljeni epsiloni posameznega elementa: eps(1)...eps(N)
   	eps_elA(:,:,inkrement+1)=eps_elA(:,:,inkrement+1)+deltae_el;

		% popravljene kappe posameznega elementa: kappa(1)...kappa(M)
   	kappa_elA(:,:,inkrement+1)=kappa_elA(:,:,inkrement+1)+deltak_el;
      
      
      for element=1:st_elementov

         norma(element)=max([norm(deltau_el(:,element),2),...\n
  		      norm(deltaR_el(:,element),2),...\n
        		norm(deltae_el(:,element),2),...\n
	         norm(deltak_el(:,element),2)]);
   	end %element=1:st_elementov
   
      disp(sprintf(' -> %2i %16.4e %15.4e %15.4e %15.3e',...\n
   	   iter,norm(deltaRu_konst,2),norm(f_konst,2),max(norma),detK_konst(inkrement+1)));

      
      % ***********************************************************
      %  prekinitev iteracije, ko je dosežena zahtevana natanènost
      % ***********************************************************
      if max([norm(deltaRu_konst,2),norm(f_konst,2)]) < 5e-8 & divergenca==0
         
         % determinanta konstrukcije    
			if ((detK_konst(inkrement) > lambda_mej) & (detK_konst(inkrement+1) < lambda_mej)) | bisekcija==1   

            bisekcija=1;
            
            stevb=stevb+1;
            
            % izpolnjena natanènost bisekcije glede determinante
            if ink_cas1 < odmik & detK_konst(inkrement+1) > lambda_mej
               
               % determinanta je enaka zahtevani vrednosti
               detK_konst(inkrement+1)=lambda_mej;
               
               % preèni prerez je spremenil definitnost
               %defc_pr(elem,sto,inkrement+1)=0;
               
               
               bisekcija=0;
               
               % ********************************************************
					%  Vrednosti kolièin v integracijskih toèkah po elementih
				   % ********************************************************

               % popravljene vrednosti specifiène spremembe dolžine in psevdoukrivljenosti
   				eps_stari=eps_ena;
   				kappa_stari=kappa_ena;
	
				   % popravljene vrednosti spremenljivke alfa za leseni prerez		
					alfaw1=alfa_w1;

				   % popravljene vrednosti meje plastiènosti v integracijskih toèkah lesenega prereza
				   sigmawy1=sigma_wy1;
				   % popravljene vrednosti elastiènega modula lesa
				   EtT=Et_T;
                   EcT=Ec_T;
             
               break
            else % ne dovolj velika natanènost bisekcije
               
               % determinanta konstrukcije
               if (detK_konst(inkrement+1) > lambda_mej) % nov izhodišèni inkrement
               
						% zaèetni posplošeni vozlišèni pomiki konstrukcije 
                  Ru_konstGKSAp=Ru_konstGKSA(:,inkrement+1);
                  
                  % posplošeni vozlišèni pomiki posameznega elementa
      				u_elAp=u_elA(:,:,inkrement+1);

						% posplošeni notranji pomiki posameznega elementa: N(0),Q(0),M(0)
						R_elAp=R_elA(:,:,inkrement+1);  

						% epsiloni posameznega elementa: eps(1)...eps(N)
				   	eps_elAp=eps_elA(:,:,inkrement+1);

						% kappe posameznega elementa: kappa(1)...kappa(M)
                  kappa_elAp=kappa_elA(:,:,inkrement+1);
                  
                  casAp=casA(inkrement);
                  
                  
                  sigma_w1p=sigma_w1(:,:,:,:,:,:,inkrement+1);
                  qu_w1p=qu_w1(:,:,:,:,:,:,inkrement+1);
                  % popravljene vrednosti mehanskih deformacij v integracijskih toèkah lesenega profila
                  D_w1p=D_w1(:,:,:,:,:,:,inkrement+1);
                  D_wp1p=D_wp1(:,:,:,:,:,:,inkrement+1);
	            	Ds_w1p=Ds_w1(:,:,:,:,:,:,inkrement+1);
                    Dc_w1p=Dc_w1(:,:,:,:,:,:,inkrement+1);
                    Dms_w1p=Dms_w1(:,:,:,:,:,:,inkrement+1);
                
                  Dt_wp=Dt_w(:,:,:,:,:,inkrement+1);
                  
  
	   				% popravljene vrednosti specifiène spremembe dolžine in psevdoukrivljenosti
   					eps_stari=eps_ena;
	   				kappa_stari=kappa_ena;
	

				   	% popravljene vrednosti spremenljivke alfa za leseni prerez		
					alfaw1=alfa_w1;

				   % popravljene vrednosti meje plastiènosti v integracijskih toèkah lesenega prereza
				   sigmawy1=sigma_wy1;
				   % popravljene vrednosti elastiènega modula lesa
				   EtT=Et_T;
                   EcT=Ec_T;

               else %(detK_konst(inkrement+1) > lambda_mej)
                  
                  if stevb == 1
                     disp(sprintf('\n -> Bisekcija do izbranega obtežnega faktorja!')); 

							% zaèetni posplošeni vozlišèni pomiki konstrukcije 
                     Ru_konstGKSAp=Ru_konstGKSA(:,inkrement);
                     
                     % posplošeni vozlišèni pomiki posameznega elementa
	      				u_elAp=u_elA(:,:,inkrement);

							% posplošeni notranji pomiki posameznega elementa: N(0),Q(0),M(0)
							R_elAp=R_elA(:,:,inkrement);  

							% epsiloni posameznega elementa: eps(1)...eps(N)
					   	eps_elAp=eps_elA(:,:,inkrement);

							% kappe posameznega elementa: kappa(1)...kappa(M)
                     kappa_elAp=kappa_elA(:,:,inkrement);
                     
                     casAp=casA(inkrement-1);
                     
                  end %stevb == 1  
                  
               end   %(detK_konst(inkrement+1) > lambda_mej)
               
               ink_cas1=ink_cas1/2;
               disp(sprintf('\n   %i. popravljeni inkrement èasa = %.3e (bisekcija)',stevb,ink_cas1))
               break
            end  % ne dovolj velika natanènost bisekcije 

			else  %(detK_konst(inkrement) < lambda_mej) &
            %(detK_konst(inkrement+1) > lambda_mej) | bisekcija==1
   
  				% popravljene vrednosti specifiène spremembe dolžine in psevdoukrivljenosti
   			eps_stari=eps_ena;
	   		kappa_stari=kappa_ena;

		   	% popravljene vrednosti spremenljivke alfa za leseni prerez		
			alfaw1=alfa_w1;

				   % popravljene vrednosti meje plastiènosti v integracijskih toèkah lesenega prereza
				   sigmawy1=sigma_wy1;
				   % popravljene vrednosti elastiènega modula lesa
				   EtT=Et_T;
                   EcT=Ec_T;

            break % prekine se zanka iter=1:15
         end 
         
      else 
         if max([norm(deltaRu_konst,2),norm(f_konst,2)]) < 0.1
            zacni=1; % sestavjati se priène celotna togostna matrika konstrukcije (raèun detK)
         end
         %prekinitev iteracije, ko je konvergenca zelo poèasna (veè kot 8 korakov)
         if max([norm(deltaRu_konst,2),norm(f_konst,2)]) < 0.1
            iter_max=10;
         else
            iter_max=8;
         end   
         if inkrement > 1 & iter > iter_max  | divergenca==1
         	divergenca=1;
            % velikost inkrementa loène dolžine zmanjša na polovico  
            ink_cas1=ink_cas1/4;%ink_cas1/5;
            stev=stev+1;
        	   disp(sprintf('\n   %i. popravljeni inkrement èasa = %.3e  (inkr.: %i)',stev,ink_cas1,inkrement))
     	      if ink_cas1 < ink_cas/(4^10)%(5^8)
  	            error(' ----- "ink_cas" je premajhen! -----')
            end
           	break
        	end   
		end
      
   end %iter=1:25
   
   

   
   % èe je divergenca=0 nadaljujemo z raèunom sicer razpolovimo ds
   if divergenca==0 & bisekcija==0

      if inkrement==(zacetek+st_korakov)
	   % izpis posplošenih vozlišènih pomikov konstrukcije
		disp(sprintf('\n **************** VOZLIŠÈNI PREMIKI KONSRTUKCIJE (cm) ****************'))
		disp(sprintf('\n%10s %13s %18s %17s','vozlišèe','pomik X','pomik Z','zasuk Y'))
		disp(' ====================================================================')
		for i=1:st_vozlisc
	   	disp(sprintf('%6i %18.7f %18.7f %19.11f',i,Ru_konstGKSA((3*(st_elementov+i-1)+1):3*(st_elementov+i),inkrement+1)'))
		   disp(' --------------------------------------------------------------------')
	   end
   end

   
   
   fprintf(fid,'\n\n  --------------------------------------------------------------');
   fprintf(fid,'\n   ->  Raèunski korak = %i  (t=%.2f h)',inkrement,casA(inkrement)/3600);
   fprintf(fid,'\n  --------------------------------------------------------------\n');
      
   
   disp(sprintf('\n   Precni prerezi z negativno definitnostjo:'))
   disp('  -------------------------------------------')
   fprintf(fid,'\n\n    Precni prerezi z negativno definitnostjo:');
   fprintf(fid,'\n   -------------------------------------------');
   st=0;
   for m=1:st_elementov
      for n=1:stopnja
         if defc_pr(m,n,inkrement+1)==-1
            disp(sprintf('    %3i. element, x = %.2f',m,dolzina(m)/2*(xg(n)+1)))
            fprintf(fid,'\n     %3i. element, x = %.2f',m,dolzina(m)/2*(xg(n)+1));
            st=st+1;
         end
      end
   end
   if st==0
      disp(sprintf('                   /'))
      fprintf(fid,'\n                    /');
   end
      
   
   
   % izpis premikov konstrukcije na datoteko IZPIS
	fprintf(fid,'\n\n\n                 *************************************\n');
	fprintf(fid,'                  VOZLISCNI PREMIKI KONSTRUKCIJE (cm)\n');  
	fprintf(fid,'                 *************************************\n');
	fprintf(fid,'\n%10s %13s %18s %17s\n','vozlisce','pomik X','pomik Z','zasuk Y');
	fprintf(fid,' ====================================================================\n');
	for i=1:st_vozlisc
	   fprintf(fid,'%6i %18.7f %18.7f %19.11f\n',i,Ru_konstGKSA((3*(st_elementov+i-1)+1):3*(st_elementov+i),inkrement+1)');
	   if i~=st_vozlisc
	      fprintf(fid,' --------------------------------------------------------------------\n');
	   end   
	end
   fprintf(fid,' ====================================================================\n');
   
   
   %izpis posplošenih notranjih pomikov elementov
   fprintf(fid,'\n\n                 *************************************\n');
	fprintf(fid,'                  POSPLOŠENI NOTRANJI POMIKI (kN,cm)\n');  
	fprintf(fid,'                 *************************************\n');
	fprintf(fid,'\n%9s %13s %16s %18s\n','element','R1(x=0)','R2(x=0)','M(x=0)');
	fprintf(fid,' ====================================================================\n');
	for i=1:st_elementov
      fprintf(fid,'%6i %16.5f %16.5f %19.5f\n',i,R_elA(:,i,inkrement+1));
      if i~=st_elementov
         fprintf(fid,' --------------------------------------------------------------------\n');
      end   
   end
   fprintf(fid,' ====================================================================\n');
   
   
   
   % ************************************************************
	%  Raèun premikov in notr. statiènih kolièin vzdolž elementov 
   % ************************************************************
   %if inkrement==(zacetek+st_korakov)
      
		for element=1:st_elementov
			   
			% preèni prerez elementa
			if isempty(pr)
		   	prerez=1;
			else   
	  			prerez=pr(element,2);
	      end
	   
   	   % ************************************************************
			%  Raèun premikov in notr. statiènih kolièin vzdolž elementov 
         % ************************************************************  
         if racunanje==1
            
         end   
                  
         if racunanje==2
            % element s konstantnim potekom epsilona in kappe
            if isempty(intersect(element,kap_const))==0
               
            else % interpolacija epsilona in kappe z Lagrangevimi polinomi  
		   	    [koord,eps(element,:),kap(element,:),u(element,:),w(element,:),fi(element,:),...\n
   		   	    silaN(element,:),silaQ(element,:),momentM(element,:)]=...\n
      		   	    NSK_gmnles(element,st_prerezov,prerez,dolzina(element),polinom,stopnja,xg,wg,...\n
         	        prirast(:,element)',G_modul(prerez),G_ploscina(prerez),...\n
            	    obt_faktor1(inkrement+1),obtn_elem(element,:),obtsp_elem(element,:),...\n
                    st_odsekov);
            end   
        end   
  
       

	      % ************************************************************
			%  Izpis premikov in notr. statiènih kolièin vzdolž elementov 
      	% ************************************************************
      	if element==1
	      	fprintf(fid,'\n\n        ********************************************************\n');
				fprintf(fid,'         PREMIKI IN NOTRANJE STATICNE KOLICINE VZDOLZ ELEMENTOV \n');  
   		   fprintf(fid,'        ********************************************************\n');
         end
         
     		obt=1; 
         % geometrijske in deformacijske kolièine konènega elementa
	      deform1(inkrement+1,(st_odsekov+1)*(element-1)+1:(st_odsekov+1)*element)=koord./100;
         deform2(inkrement+1,(st_odsekov+1)*(element-1)+1:(st_odsekov+1)*element)=...\n
            u(element,:);
         deform3(inkrement+1,(st_odsekov+1)*(element-1)+1:(st_odsekov+1)*element)=...\n
            w(element,:); 
         deform4(inkrement+1,(st_odsekov+1)*(element-1)+1:(st_odsekov+1)*element)=...\n
            1000*eps(element,:); % v promilih
         deform5(inkrement+1,(st_odsekov+1)*(element-1)+1:(st_odsekov+1)*element)=...\n
            1000*kap(element,:); % v promilih
         
         
         
         % risanje upogibnih momentov konstrukcije (mater. in geom. linearnost)
	      momenti(inkrement,(st_odsekov+1)*(element-1)+1:(st_odsekov+1)*element)=...\n
            momentM(element,:)./100;
               
      	% risanje osnih sil konstrukcije (mater. in geom. linearnost)
         osnesile(inkrement,(st_odsekov+1)*(element-1)+1:(st_odsekov+1)*element)=...n
            silaN(element,:);
      
      	% risanje preènih sil konstrukcije (mater. in geom. linearnost)
         precnesile(inkrement,(st_odsekov+1)*(element-1)+1:(st_odsekov+1)*element)=...\n
            silaQ(element,:);
         
   
         
         % izpis premikov in notranje statiènih kolièin vzdolž elementov na datoteko IZPIS
			fprintf(fid,'\n\n         ****************************** %i. ELEMENT ******************************\n',element);
         
         fprintf(fid,'\n%60s','D(prom.):');
         fprintf(fid,'\n%7s %13s %13s %17s','x(m)','epsilon','kappa','sp.rob');
         fprintf(fid,'%10s','zg.rob');
         fprintf(fid,'\n ===================================================================================================\n');
         for st=1:st_odsekov+1
            fprintf(fid,'%7.2f %16e %15e %12.3f',[koord(1,st)./100 eps(element,st) kap(element,st) 1000*(eps(element,st)+kap(element,st)*(sum(visina(prerez,:))-zT(prerez)))]);	      
            fprintf(fid,'%10.3f\n',1000*(eps(element,st)+kap(element,st)*(-zT(prerez))));
         end   
         fprintf(fid,' ===================================================================================================\n');
      
			fprintf(fid,'\n%7s %10s %10s %10s %11s %11s %13s','x(m)','u(cm)','w(cm)','fi()','N(kN)','Q(kN)','M(KNcm)');
			fprintf(fid,'\n =======================================================================================\n');
	   
		   fprintf(fid,'%7.2f %10.4f %10.4f %11.5f %10.2f %10.2f %14.2f\n',[koord'./100 u(element,:)' w(element,:)' fi(element,:)' silaN(element,:)' silaQ(element,:)' momentM(element,:)']');
         fprintf(fid,' =======================================================================================\n\n');
         

		end % element=1:st_elementov
	%end%inkrement=(zacetek+st_korakov)
   
      
   end % divergenca==0 & bisekcija==0
  
   
   end % while divergenca==1
   
   
end %for inkrement=(zacetek+1):(zacetek+st_korakov)

fclose(fid);


% porabljen cas za racun konstrukcije
disp(sprintf('\n  Porabljen èas: %3.2f sek',cputime-tm))


   


