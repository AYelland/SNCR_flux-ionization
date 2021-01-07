% GALACTIC COSMIC RAY FLUX

  clear variables
  clear all
  close all
  format long
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-- Background GCR Flux --%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% This section of the code calculated the Background Galactic Cosmic Ray (GCR) Spectrum

  GCRfile = csvread("SN_50pc_CaseB.csv") ;
    tmaxGCR = size(GCRfile, 2)-1 ;
    EmaxGCR = size(GCRfile, 1)-1 ;
    
  % Time Spectrum
    for j = 1:1:tmaxGCR
      tdatGCR(j,1) = GCRfile(1,j+1) ;
    endfor

  % Energy Spectrum
    for i = 1:1:EmaxGCR
      EdatGCR(i,1) = GCRfile(i+1,1) ;
    endfor

  % GCR Flux Data from Input File
    for i = 1:1:EmaxGCR
      for j = 1:1:tmaxGCR
        fluxdatGCR(i,j) = GCRfile(i+1,j+1) ; % Not used in calculations
      endfor
    endfor

  % GCR Flux Background Spectrum
    % Equation from Haungs: Rebel and Roth: 2003: Reports on Progress in Physics: vol. 66: p. 1145-1206
      % Two power laws: one below the "knee" & one above the "knee" (knee location = 3.2E6 GeV)
      % Data for formula, page 1179   
      
    fluxGCRbackground = zeros(EmaxGCR,1) ;

    f2 = -2.68 ;     %(below the "knee")    >>> ( −(2.68 ± 0.06) )
    f3 = -3.06 ;     %(above the "knee")    >>> ( −(3.06 ± 0.08) )
    Eknee = 3.2E6 ;  %(at the "knee")       >>> (knee = 3.2E15 ± 1.2E15 eV = 3.2E6 ± 1.2E6 GeV)
    
    f1_below = 1.8E4 %from Dr. Andrew Overholt's email 15March2016
    f1_above = f1_below*(Eknee^(f2-f3)) %computed from formula for straight line in logspace [ y = b*(x^m) ]
    
    for ecur = 1:1:EmaxGCR
    %printf("EdatGCR(%i) = %f \n", ecur , EdatGCR(ecur))
      if (EdatGCR(ecur) < Eknee)
      %printf("EdatGCR(%i) = %f \n", ecur , EdatGCR(ecur))
        fluxGCRbackground(ecur) = f1_below*(EdatGCR(ecur))^f2 ;
        %disp("Below Knee")
        %printf("EdatGCR(%i) = %f , fluxGCRbackground(%i) = %f \n", ecur , EdatGCR(ecur), ecur, fluxGCRbackground(ecur)*EdatGCR(ecur)*EdatGCR(ecur))
      else
        fluxGCRbackground(ecur) = f1_above*(EdatGCR(ecur))^f3 ;
        %disp("Above Knee")
        %printf("EdatGCR(%i) = %f , fluxGCRbackground(%i) = %f \n", ecur , EdatGCR(ecur), ecur, fluxGCRbackground(ecur)*EdatGCR(ecur)*EdatGCR(ecur))
      endif
      fluxGCRbackgroundPlot(ecur) = fluxGCRbackground(ecur)*EdatGCR(ecur)*EdatGCR(ecur) ; % Used for plotting purposes
      %disp(" ")
    endfor
    
    
    % GCR Background Flux Plotting

      figure(1)
      
      loglog((EdatGCR),(fluxGCRbackgroundPlot), "k", "linewidth", 2)
        title("Background GCR Flux")
        xlabel("Proton Energy [GeV]")
        ylabel("CRFlux * E^2 [(protons/(GeV m^2 s sr))*(GeV)^2]")
##          axis([10^6,10^7,10^4,10^6])
        grid on

##           figfile = "Plot - GCR_Flux" ;
##           print(figfile,"-dpng") ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-- SNCR Flux --%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Read-in CR Flux Data (This is the "variable" data we compare with.)

  % This data represents the amount of cosmic rays that enter the atmosphere due to the specified supernova 
    % Utilizes same "tdat" & "Edat" values as from GCR Flux
    % (Energy Bins "Edat" vs. Time "tdat")

    SNCRfile = csvread("SNCR_flux_50pc_Melott_CaseB.csv") ; ;
      tmax = size(SNCRfile, 2)-1 ;
      Emax = size(SNCRfile, 1)-1 ;
      
    % Time Spectrum
      for j = 1:1:tmax
        tdat(j,1) = SNCRfile(1,j+1) ;
      endfor
    
    % Energy Spectrum
      for i = 1:1:Emax
        Edat(i,1) = SNCRfile(i+1,1) ;
      endfor

    % SNCR Flux Data from Input File
      for i = 1:1:Emax
        for j = 1:1:tmax
          fluxdat(i,j) = SNCRfile(i+1,j+1) ;
          fluxdatPlot(i,j) = fluxdat(i,j)*((Edat(i))^2) ;
        endfor
      endfor   

  % Flux Plotting
        
        labels = cellstr(["100yr";"300yr";"1000yr";"3000yr";"10000yr";"30000yr";"100000yr";"300000yr";"1000000yr";"3000000yr";"GCR"]) ;
        figure(2)
        
        hold on
        
        loglog(Edat, fluxdatPlot, "linewidth", 2)
        loglog(EdatGCR,fluxGCRbackgroundPlot, "k", "linewidth", 2) 
 
          legend(labels, "location", "northeastoutside")
          title("Supernova Cosmic Ray Flux")
          axis([10,10^(7),10^(-6),10^(6)])            
          xlabel("Energy [GeV]")
          ylabel("CR Flux [(protons/(GeV m^2 s sr))*(GeV^2)]")
          grid on
        
        hold off


    
    

