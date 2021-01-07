 
% SUPERNOVA COSMIC RAY ATMOSPHERIC FLUX & IONIZATION 
% (50 parsecs away -- set your distance below in "Supernova CR flux" section)
%
% For background info, see publications:
   % Thomas et. al, ApJ Letters, 826, 1, 2016 
   %  https://arxiv.org/abs/1605.04926
   % Melott et al., ApJ, 840, 105, 2017
   %  https://arxiv.org/abs/1702.04365
%
% Created by Alex Yelland and Brian Thomas, Washburn University, 2020
   % Shared under Creative Commons license:
   %  Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
   %  https://creativecommons.org/licenses/by-nc-sa/4.0/
%
% Questions should be directed to Brian Thomas, brian.thomas@washburn.edu

  clear variables
  clear all
  close all
  format long

% -------------------- MAIN EDITING SECTION ----------------------------


  % Input/output files
  % When changing parameters (such as distance from supernova), change filenames accordingly

    % GCR Flux Files:
      gcrfluxfile_input = "SN_50pc_CaseB.csv" ;               %(Input File)

    % SNCR Flux Files:
      sncrfluxfile = "SNCR50pc-flux.csv" ;                    %(Output File)

    % Ionization Data Files:
      sncrfluxReadin = "SNCR50pc-flux.csv" ;                  %(Output File)
      gcrfluxfile_output = "GCR-flux.dat" ;                   %(Output File)
      gcrfluxplotfile = "GCR-flux_Plot.dat" ;                 %(Input File)
      Dionsfile = "Dimitra_Ionization.csv" ;                  %(Input File)

    % Atmospheric Ionization Files:
      atmofile = "AtmoDepth&Density_GroundUp.csv" ;           %(Input File)
      
      % Plot Title
        IonsPlotTitle = "Atmospheric Ionization (50pc)" ;
      
    % Ionization Data Files Out:
      sncrionsCSV = "SNCR50pc-ions_cm3.csv" ;                 %(Output File)
      sncrionsDAT = "SNCR50pc-ions_cm3.dat" ;                 %(Output File)
      sncrionsDAT_GCR_cm2 = "GCR-ions_cm2.dat" ;              %(Output File)
      sncrionsDAT_GCR_cm3 = "GCR-ions_cm3.dat" ;              %(Output File)
      sncrionsDAT_SNCR_cm2 = "SNCR50pc-ions_cm2-%i.dat" ;     %(Output File)
      sncrionsDAT_SNCR_cm3 = "SNCR50pc-ions_cm3-%i.dat" ;     %(Output File)
      
      

  % Logical Code Variables to Determine which Sub-Codes to Run { "true/run" (1) or "false/skip" (0) }

    flux                    = 1 ;
      fluxGCRPlot           = 0 ;
      fluxSNCRPlot          = 0 ;
      fluxOutput            = 1 ;
      
    ions                    = 1 ;
      readin_fluxSNCRPlot   = 0 ; %valid only when flux = 0
      ionsDimitraPlot       = 1 ;
      ionsInterpPlots       = 0 ;
      ionsPlot              = 1 ;
      ionsOutput            = 1 ;

  % Plotting Values/Properites

    labels = {"GCR";"100yr";"300yr";"1000yr";"3000yr";"10000yr";"30000yr";"100000yr";"300000yr";"1000000yr";"3000000yr"} ;
    colors = [ 0.0,0.0,0.0 ; 0.9,0.0,0.0 ; 0.9,0.5,0.2 ; 1.0,0.8,0.0 ; 0.3,0.6,0.0 ; 0.0,1.0,0.0 ; 0.3,0.4,1.0 ; 0.5,0.8,1.0 ; 0.6,0.2,0.9 ; 1.0,0.4,1.0 ; 0.6,0.3,0.0 ] ;


% -------------------------- END MAIN EDITING SECTION ------------------------- 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-- Cosmic Ray Flux --%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if(flux == true)
  disp("---------- FLUX -----------")


  % -------------------------- Background GCR Flux ----------------------------    
    %% This section of the code calculated the Background Galactic Cosmic Ray (GCR) Spectrum
    disp("GCR Flux:: ")
    
    cd("ReadinFiles") 
    GCRfile = csvread(gcrfluxfile_input) ;
    cd ..
      tmaxGCR = size(GCRfile, 2)-1 ; %yrs
      EmaxGCR = size(GCRfile, 1)-1 ; %GeV
      
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

      f1_below = 1.8E4 ; %(m2*sr*s)^(-1) %From Tanabashi et al. (Particle Data Group), Phys. Rev. D 98, 030001 (2018), DOI: 10.1103/PhysRevD.98.030001, Chapter 29 (Cosmic Rays)
      f1_above = f1_below*(Eknee^(f2-f3)) ; %computed from formula for straight line in logspace [ y = b*(x^m) ]

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
        fluxGCRbackgroundPlot(ecur) = fluxGCRbackground(ecur)*EdatGCR(ecur)*EdatGCR(ecur) ; %Used for plotting purposes

        GCRfluxFileOutput(ecur,1) = fluxGCRbackground(ecur) ; %(m2*sr*s*GeV)^(-1)
        GCRfluxPlotFileOutput(ecur,1) = EdatGCR(ecur) ; %GeV
        GCRfluxPlotFileOutput(ecur,2) = fluxGCRbackgroundPlot(ecur) ; %(m2*sr*s*GeV)^(-1) * (GeV)^2
        %disp(" ")
      endfor
      
      
      % GCR Background Flux Plotting
      
        if((fluxGCRPlot) == true)
        disp("  GCR Flux Plot: ")

        figure(1)
          
          loglog((EdatGCR),(fluxGCRbackgroundPlot), "k", "linewidth", 1)
            title("Background GCR Flux")
            xlabel("Proton Energy [GeV]")
            ylabel("CRFlux * E^2 [(protons/(GeV m^2 s sr))*(GeV)^2]")
            legend("GCR", "location", "northeastoutside")
            grid on

          cd("Plots")
          figfile = "Plot - GCR_Flux" ;
          print(figfile,"-dpng") ;
          cd ..

        figure(2)
          
          loglog((EdatGCR*((10)^9)),(fluxGCRbackground), "k", "linewidth", 1)
            title("Background GCR Flux")
            xlabel("Proton Energy [eV]")
            ylabel("CRFlux [protons/(GeV m^2 s sr)]")
            legend("GCR", "location", "northeastoutside")
            axis([10^(8),10^(21),10^(-29),10^(4)])
            grid on

          cd("Plots")
          figfile = "Plot - GCR_Flux2" ;
          print(figfile,"-dpng") ;
          cd ..

        endif %plotting


  % -------------------------- Supernova CR Flux ----------------------------    
    % This section of the code is used to develop the cosmic ray spectrum modified from Melott2017 data (2020Feb04 AY)
    disp("SNCR Flux: ")


    %Constants

      % Speed of Light
          c    = 2.99792458e+8 ;                                  % m/s
          c_cm = c*100 ;                                          % cm/s
          c_pc = (c*(365*24*60*60))/(30856775812799585) ;         % parsec/year

      % Initial Diffusion Rate
          Dnot = 2.0e+28 ;                                        % cm^(2)/s

      % Initial Energy
          Enot = 1.0 ;                                            % GeV

      % Cutoff Energy
          Ec = 1.0e+6 ;                                           % GeV    

      % Source Distance
          r_pc = 50 ;                                             % parsec (pc)
          r_cm = r_pc*(30856775812799585*100) ;                   % cm 
      
    % Time Spectrum

        tdat_yrs = tdatGCR ;
        tdat = (tdat_yrs + (r_pc/c_pc))*(365*24*60*60) ;          % s
          % Note on Time:  from MK email 11Sept2019
            % tdat = 100 yr means tdat = 100 yr + r/c 
            % Where r = distance in parsecs and c = speed of light in parsecs / year
            % So for the "100 year" calculations, the t value in the equation should actually be about 267 years (but converted to seconds).
        tmax = size(tdat, 1) ;

    % Energy Spectrum

        Edat = EdatGCR ;                                          % GeV 
        Emax = size(Edat, 1) ;

    % Computing the CR Flux/Intensity

      for i = 1:1:tmax
        for j = 1:1:Emax

        % Source Spectrum
          factor = 3.0e+52 ;      %from MK email 8Sept2019 
          Q(j) = factor*(Edat(j).^(-2.2)).*exp(-Edat(j)/Ec) ;     % protons/GeV 

        % Diffusion Approximation
          D(j) = (Dnot).*(Edat(j)./Enot).^(1/3) ;                 % cm^(2)/s

        % Effective Diffusion Distance
          rdiff(j,i) = sqrt(4.*D(j).*tdat(i)) ;                   % cm

        % Melott SN CR Density Equation
          n(j,i) = (Q(j)./((pi^(3/2)).*(rdiff(j,i).^(3))).*(exp(-r_cm^(2)./rdiff(j,i).^(2)))) ;          % protons/((GeV)(cm^(3)))

        % CR Intensity/Flux
          fluxdat_cm2(j,i) = (c_cm/(4*pi))*(n(j,i)) ;            % (protons/((GeV)(cm^(2))(s)(sr)))
          fluxdat_m2(j,i) = fluxdat_cm2(j,i)*100*100 ;           % (protons/((GeV)(m^(2))(s)(sr)))
          fluxdatPlot(j,i) = fluxdat_m2(j,i)*((Edat(j))^2) ;     % (protons/((GeV)(m^(2))(s)(sr)))*(GeV^2)

        endfor
      endfor


    % SNCR Flux Plotting (w/ GCR Background Flux)
    
      if((fluxSNCRPlot) == true)
      disp("  SNCR Flux Plot: ")
                
        figure(3)
        
          hold on
          loglog(EdatGCR, fluxGCRbackgroundPlot, "color", colors(1,:), "linewidth", 1)
          for i = 1:1:tmax
            loglog(Edat, fluxdatPlot(:,i), "color", colors(i+1,:), "linewidth", 1)
          endfor
          hold off
          
          legend(labels, "location", "northeastoutside")
          title("Supernova Cosmic Ray Flux")
          axis([10,10^(7),10^(-6),10^(6)])            
          xlabel("Energy [GeV]")
          ylabel("CR Flux [(protons/(GeV m^2 s sr))*(GeV^2)]")
          grid on
        
        cd("Plots")    
        figfile = "Plot - SNCR_Flux" ;
        print(figfile,"-dpng") ;
        cd ..

      endif %plotting


  % -------------------------- Flux Output ----------------------------    
      
    if ((fluxOutput) == true)
    disp("Output Flux Data: ")
    
      cd("OutputFiles")
      
      % SNCR Flux
        SNCRfluxFileOutput = zeros(Emax+1,tmax+1) ;
        for i = 1:1:Emax
          for j = 1:1:tmax
            SNCRfluxFileOutput(i+1,1) = Edat(i) ; %GeV
            SNCRfluxFileOutput(1,j+1) = tdat_yrs(j) ; %yrs
            SNCRfluxFileOutput(i+1,j+1) = fluxdat_m2(i,j) ; %(protons/((GeV)(m^(2))(s)(sr)))
          endfor
        endfor
      
        dlmwrite(sncrfluxfile, SNCRfluxFileOutput, ",") ;
      
      % GCR Flux
        save ("-ascii", gcrfluxfile_output, "GCRfluxFileOutput") ;
        save ("-ascii", gcrfluxplotfile, "GCRfluxPlotFileOutput") ;

      cd ..
    
    endif %fluxOutput

  endif %flux


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-- Ions --%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% This section of the code is used to determine the atmospheric ionization caused by the cosmic ray flux

  if(ions == true)
  disp("---------- IONS -----------")
  

  % -------------------------- Read-in CR Flux Data ----------------------------
  if((flux) != true)
  disp("SNCR Flux Read-In:")  

    % Read-in CR Flux Data
    % This data represents the amount of cosmic rays that enter the atmosphere due to the specified supernova & background GCR
    
      cd("ReadinFiles") 
      SNCRfile = csvread(sncrfluxReadin) ;
      cd ..
      % (Energy Bins "Edat" vs. Time "tdat_yrs")
        tmax = size(SNCRfile, 2)-1 ;
        Emax = size(SNCRfile, 1)-1 ;
      
      cd("ReadinFiles") 
      GCRFile = load(gcrfluxplotfile) ; %computed and created above
      cd ..
      % (column 1 = Energy Bins "EdatGCR" , column 2 = GCR Flux Data * E^2 "fluxGCRbackgroundPlot")
        EmaxGCR = size(GCRFile, 1) ;
        
      % Time Spectrum
        for j = 1:1:tmax
          tdat_yrs(j,1) = SNCRfile(1,j+1) ;
        endfor
      
      % Energy Spectrums
        for i = 1:1:Emax
          Edat(i,1) = SNCRfile(i+1,1) ;
        endfor
        
        for i = 1:1:EmaxGCR
          EdatGCR(i,1) = GCRFile(i,1) ;
        endfor

      % Flux Data from Input Files
        for i = 1:1:Emax
          for j = 1:1:tmax
            fluxdat_m2(i,j) = SNCRfile(i+1,j+1) ; % (protons/((GeV)(m^(2))(s)(sr)))
            fluxdatPlot(i,j) = fluxdat_m2(i,j)*((Edat(i))^2) ; % (protons/((GeV)(m^(2))(s)(sr)))*(GeV^2)
          endfor
        endfor
        
        for i = 1:1:EmaxGCR
          fluxGCRbackgroundPlot(i,1) = GCRFile(i,2) ;
        endfor
     
    
    % Flux Plotting
    
        if((readin_fluxSNCRPlot) == true)
        disp("  SNCR Flux Read-In Plot: ")
                    
          figure(4)
         
            hold on
            loglog(EdatGCR,fluxGCRbackgroundPlot, "color", colors(1,:), "linewidth", 1)
            for i = 1:1:tmax
              loglog(Edat, fluxdatPlot(:,i), "color", colors(i+1,:), "linewidth", 1)
            endfor
            hold off
            
            legend(labels, "location", "northeastoutside")
            title("Supernova Cosmic Ray Flux - OLD")
##            axis([10^(6),10^(9),10^(-55),10^(5)])
            axis([10,10^(7),10^(-6),10^(6)])            
            xlabel("Energy [GeV]")
            ylabel("CR Flux [(protons/(GeV m^2 s sr))*(GeV^2)]")
              set(gca, "ytick", 10.^(-55:5:5))
            grid on
          
          cd("Plots")   
          figfile = "Plot - SNCR_Flux (old)" ;
          print(figfile,"-dpng") ;
          cd ..

        endif %readin_fluxSNCRPlot
            
  endif %fluxSNCR read-in
    

    
  % ------------------- Read-in Dimitra's Ionization Data ----------------------
  disp("DIons Read-In:")  

    % Read-in Dimitra's Ionization Data: { D = Dimitra }
      % This data represents the ionization of a single cosmic ray (a proton) entering the atmosphere causing the ionizing "shower" effect.
      % (Atmospheric Pressure Bins "PDat" vs. Energy Bins "DEdat") 

      cd("ReadinFiles") 
      dimitraDat = csvread(Dionsfile) ;
      cd ..
        % Original file titled "merged-table_Dimitra-Ions_march2016.csv"
        % Downloaded from ftp://kusmos.phsx.ku.edu/data/melott/crtables/atmoionization/ on April 21, 2020
      
        DEmax = size(dimitraDat, 2) - 1 ;     % the "minus one" is to compensate for the first column in the formating of the file
        DEdat = zeros(DEmax,1) ;
        for i = 1:1:DEmax
          DEdat(i) = dimitraDat(1,i+1) ; %GeV
          % Contained within first row of file
        endfor

        Pmax = size(dimitraDat, 1) - 1 ;    % the "minus one" is to compensate for the first row in the formating of the file
        PDat = zeros(Pmax,1) ;
        for j = 1:1:Pmax
          PDat(j) = dimitraDat(j+1,1) ; %hPa
          % Contained within first column of file
        endfor

        DIons = zeros(DEmax,Pmax) ;
        DIons_perGeV = zeros(DEmax,Pmax) ;
        for i = 1:1:DEmax
           for j = 1:1:Pmax
              DIons(i,j) = dimitraDat(j+1,i+1) ; %(ion pairs) % Transposes table from the read-in file Dionsfile into "DIons" array
              DIons_perGeV(i,j) = DIons(i,j)/DEdat(i) ;
           endfor
        endfor



    % Finding the Altitude Pressure Bin at which the Maximum Ionization Rate Occurs (for each primary proton energy value)
    % NOTE!! The "PDat" values here (taken from the table input file) are actually *pressure* values (in hPa), not altitudes

        maxIon = zeros(DEmax,1) ;
        maxIonE = zeros(DEmax,1) ;
        maxIonP = zeros(DEmax,1) ;

        for i = 1:1:DEmax
           for j = 1:1:Pmax
              if (DIons(i,j) > maxIon(i))
                maxIon(i) = DIons(i,j) ;
                %maxIon(i) = max(DIons(i,:)) ;
                maxIonE(i) = DEdat(i) ;
                maxIonP(i) = PDat(j) ;
              endif
           endfor
        endfor

    % Dimitra Ionization Plotting

      if((ions)&(ionsDimitraPlot) == true)
      disp("  DIons Plot: ")

        figure(5)
        
        loglog(DIons,PDat, "linewidth", 1)
          title("Dimitra's Ionization")
          xlabel("Ionization [(ion pairs)/proton]")
          ylabel("Pressure [hPa]")
          %axis([,,,])
          grid on

        cd("Plots")
        figfile = "Plot - Dimitra_Ions" ;
        print(figfile,"-dpng") ;
        cd ..

      endif %plotting



  % ----------- Interpolating Energy Bins, Ionization Data, & Flux Data (GCR & SNCR) onto High-Resolution E-Grid ------------
  disp("Interpolation:")  

    InterpBinNum = 121 ;    % Use an odd number

    % Setting Higher Resolution to Energy Bins for Interpolation ("E" units: GeV) 
    disp("  Energy Bins Interp:")  
      
      % Minimum Energy Value in both Energy Ranges
        if (min(Edat) < min(DEdat))
          Ebinsmin = min(DEdat) ;
        else
          Ebinsmin = min(Edat) ;
        endif
        Elo = log10(Ebinsmin) ;
      
      % Maximum Energy Value in both Energy Ranges
        if (max(Edat) > max(DEdat))
          Ebinsmax = max(DEdat) ;
        else
          Ebinsmax = max(Edat) ;
        endif
        Ehi = log10(Ebinsmax) ;

      % How many e-grid points drop off in the beginning? (too low in E)
      
        % Edat
          rmLowE = 1 ;
          while (Edat(rmLowE) < Ebinsmin)
            rmLowE = rmLowE+1 ;
          endwhile
          rmLowE = rmLowE - 1 ; % (the minus 1 is for the indexing of the loop)
      
        % DEdat
          D_rmLowE = 1 ;
          while (DEdat(D_rmLowE) < Ebinsmin)
            D_rmLowE = D_rmLowE+1 ;
          endwhile
          D_rmLowE = D_rmLowE - 1 ; % (the minus 1 is for the indexing of the loop)
          
      % How many e-grid points drop off at the end? (too high in E)
      
        % Edat
          rmHighE = Emax ;
          while (Edat(rmHighE) > Ebinsmax)
            rmHighE = rmHighE-1 ;
          endwhile
          rmHighE = Emax - rmHighE ;

        % DEdat
          D_rmHighE = DEmax ;
          while (DEdat(D_rmHighE) > Ebinsmax)
            D_rmHighE = D_rmHighE-1 ;
          endwhile
          D_rmHighE = DEmax - D_rmHighE ;


    % Interpolation: Creating new "Ebins" grid based from new "Ebinsmin" & "Ebinsmax" boundaries

        IonsEmax = InterpBinNum ;

        Ebins = zeros(IonsEmax,1) ;        
        Ebins_Rows = logspace(Elo, Ehi, IonsEmax);
        for i = 1:1:IonsEmax
          Ebins(i,1) = Ebins_Rows(i) ; %GeV
        endfor
        % "logspace" performs the following function:
          % Esteps = IonsEmax-1 ;
          % Estepsize = (Ehi - Elo)/Esteps ;

          % for i = 1:1:IonsEmax
          %   Ebins(i) = 10.^(Elo+(Estepsize*(i-1))) ; % (Reversing the log10 used in "Elo" & "Ehi")
          % endfor
        % Ebins
         


    % Interpolation Loop: Ionization Data from "DEdat" grid to a new "Ebins" grid
    disp("  Ions Interp:")  

      newDIons = zeros(IonsEmax,Pmax) ; %(ion pairs)

      for j = 1:1:Pmax

        inew =1 ;
        iin = D_rmLowE+1 ; % index for Dimitra's "DEdat", need to start at "Ebinsmin" (10GeV) [File starts at 0.3GeV] 
        f1 = 0 ;
        f2 = 0 ;
        f3 = 0 ;

        for i = 1:1:IonsEmax

          if (Ebins(inew) == DEdat(iin))
            newDIons(inew,j) = DIons(iin,j) ;
            inew = inew + 1 ;

          elseif ((Ebins(inew) > DEdat(iin)) & (Ebins(inew) < DEdat(iin+1)))
            f1 = log10(Ebins(inew)/DEdat(iin)) ;
            f2 = log10(DEdat(iin+1)/DEdat(iin)) ;
            f3 = f1/f2 ;
            newDIons(inew,j) = (DIons(iin+1,j)^f3)*(DIons(iin,j)^(1-f3)) ;
            inew = inew +1 ;
          endif

          if (inew > IonsEmax)
            break 
          endif

          if (iin >= (DEmax-D_rmHighE))
            break
          endif

          if (Ebins(inew) >= DEdat(iin+1))
            iin = iin + 1 ;
          endif

        endfor
      endfor
      %newDIons

      newDIons_perGeV = zeros(IonsEmax,Pmax) ; %(ion pairs)/GeV
      for i = 1:1:IonsEmax
         for j = 1:1:Pmax
            newDIons_perGeV(i,j) = newDIons(i,j)/Ebins(i) ;
         endfor
      endfor
      %newDIons_perGeV



    % Interpolation Loop: Flux Data from "Edat" grid to a new "Ebins" grid
    disp("  SNCR Flux Interp:")  
              
        newfluxdat = zeros(IonsEmax,tmax) ; % (protons/((GeV)(m^(2))(s)(sr)))
                                            %drop data points that are beyond ionization table energy range

        fluxEmax = InterpBinNum ;
        % Note: fluxEmax = IonsEmax = InterpBinNum
        
        for j = 1:1:tmax 

          inew = 1 ;
          iin = rmLowE+1 ;
          f1 = 0 ;
          f2 = 0 ;
          f3 = 0 ;

          for i = 1:1:fluxEmax 
            
            if (Ebins(inew) == Edat(iin))
              newfluxdat(inew,j) = fluxdat_m2(iin,j) ;
              inew = inew + 1 ;

            elseif ((Ebins(inew) > Edat(iin)) & (Ebins(inew) < Edat(iin+1)))
              f1 = log10(Ebins(inew)/Edat(iin)) ;
              f2 = log10(Edat(iin+1)/Edat(iin)) ;
              f3 = f1/f2 ;
              newfluxdat(inew,j) = (fluxdat_m2(iin+1,j)^f3)*(fluxdat_m2(iin,j)^(1-f3)) ; %(protons/((GeV)(m^(2))(s)(sr)))
              inew = inew +1 ;
            endif

            if (inew > fluxEmax)
              break
            endif 

            if (iin >= (Emax-rmHighE))
              break
            endif

            if (Ebins(inew) >= Edat(iin+1))
              iin = iin + 1 ;
            endif

          endfor
        endfor

        for i = 1:1:fluxEmax
          for j = 1:1:tmax
            newfluxdatPlot(i,j) = newfluxdat(i,j)*((Ebins(i))^2) ; % (protons/((GeV)(m^(2))(s)(sr))) * (GeV)^2
          endfor
        endfor



    % Recalculating Background GCR Flux based from "Ebins"
    disp("  GCR Flux Recalculation:")  

      % Equation from Haungs: Rebel and Roth: 2003: Reports on Progress in Physics: vol. 66: p. 1145-1206
        % Two power laws: one below the "knee" & one above the "knee" (knee location = 3.2E6 GeV)

      newfluxGCRbackground = zeros(IonsEmax,1) ;

      f2 = -2.68 ;     %(below the "knee")    >>> ( −(2.68 ± 0.06) )
      f3 = -3.06 ;     %(above the "knee")    >>> ( −(3.06 ± 0.08) )
      Eknee = 3.2E6 ;  %(at the "knee")       >>> (knee = 3.2E15 ± 1.2E15 eV = 3.2E6 ± 1.2E6 GeV)

      f1_below = 1.8E4 ; %(m2*sr*s)^(-1) %From Tanabashi et al. (Particle Data Group), Phys. Rev. D 98, 030001 (2018), DOI: 10.1103/PhysRevD.98.030001, Chapter 29 (Cosmic Rays)
      f1_above = f1_below*(Eknee^(f2-f3)) ; %computed from formula for straight line in logspace [ y = b*(x^m) ]

      for ecur = 1:1:IonsEmax
      %printf("Ebins(%i) = %f \n", ecur , Ebins(ecur))
        if (Ebins(ecur) < Eknee)
        %printf("Ebins(%i) = %f \n", ecur , Ebins(ecur))
          newfluxGCRbackground(ecur) = f1_below*(Ebins(ecur))^f2 ;
          %disp("Below Knee")
          %printf("Ebins(%i) = %f , newfluxGCRbackgroundPlot(%i) = %f \n", ecur , Ebins(ecur), ecur, newfluxGCRbackground(ecur)*Ebins(ecur)*Ebins(ecur)
        else
          newfluxGCRbackground(ecur) = f1_above*(Ebins(ecur))^f3 ;
          %disp("Above Knee")
          %printf("Ebins(%i) = %f , newfluxGCRbackgroundPlot(%i) = %f \n", ecur , Ebins(ecur), ecur, newfluxGCRbackground(ecur)*Ebins(ecur)*Ebins(ecur)
        endif
        newfluxGCRbackgroundPlot(ecur) = newfluxGCRbackground(ecur)*Ebins(ecur)*Ebins(ecur) ; %(m2*sr*s*GeV)^(-1)*(GeV)^2  (Used for plotting purposes)
        %disp(" ")
      endfor
      newfluxGCRbackground ; %(m2*sr*s*GeV)^(-1)
      newfluxGCRbackgroundPlot ; %(m2*sr*s*GeV)^(-1)*(GeV)^2


    % Interpolation Plotting
    
      if((ions)&(ionsInterpPlots) == true)
      disp("  Interpolation Plots:")
      
        % Interpolated Flux Plot (MK's Expanded)
          
          figure(6)

            hold on
            loglog(Ebins, newfluxGCRbackgroundPlot, "color", colors(1,:), "linewidth", 1)
            for i = 1:1:tmax
              loglog(Ebins, newfluxdatPlot(:,i), "color", colors(i+1,:), "linewidth", 1)
            endfor
            hold off

            legend(labels, "location", "northeastoutside")
            title("Supernova Cosmic Ray Flux (HiRes)")
            axis([10,10^(6),10^(-6),10^(6)])          
            xlabel("Energy [GeV]")
            ylabel("CR Flux [(protons/(GeV m^2 s sr))*(GeV^2)]")
            grid on
          
          cd("Plots")    
          figfile = "Plot - SNCR_Flux_HiRes";
          print(figfile,"-dpng");
          cd ..
            
            
        % Interpolated Ions Plot (Dimitra's Expanded)
        
          figure(7)
        
          loglog(newDIons, PDat, "linewidth", 1)
            title("Dimitra's Ionization (HiRes)")
            xlabel("Ionization [(ion pairs)/proton]")
            ylabel("Pressure [hPa]")
            %axis([,,,])
            grid on

          cd("Plots")
          figfile = "Plot - Dimitra_Ions_HiRes" ;
          print(figfile,"-dpng") ;
          cd ..

      endif %plotting



  % ------------- Atmospheric Ionization due to Background GCR Flux ------------
  disp("Atmospheric Data:")

      % Atmospheric Depth & Density     [ Need to convert to per cm3 (from per cm2) ]
        
        cd("ReadinFiles") 
        atmoFileIn = csvread('AtmoDepth&Density_GroundUp.csv') ;
        cd ..
  
          % Notes on File Format:
          % - Altitude (column 1) is measured in (km)
          % - Atomspheric Depth (column 2) is measured in (g/cm^2)
          % - Atmospheric Density (column 3) is measured in (g/cm^3)
          % - Row 1/2  - Empty due to file read-in formating
          % - Data/Altitude moves from ground to sky { Dimitra's Ionization Table "newDIons" (origianlly from "DIons") is from sky to ground }

        atmoAltmax = (size(atmoFileIn,1) -2) ; % The "minus two" is for the formating 
        % NOTE: atmoAltmax = Pmax 
        
        for j = 1:1:size(atmoFileIn,2)
          for i = 1:1:atmoAltmax
            atmoConvArr(i,j) = atmoFileIn(i+2,j) ;
          endfor
        endfor
        
        for i = 1:1:atmoAltmax
          atmoAlt(i,1) = atmoConvArr(i,1) ;
          atmoDepth(i,1) = atmoConvArr(i,2) ;
          atmoDensity(i,1) = atmoConvArr(i,3) ;
        endfor
        
  disp("Ionization Calculation:")
      % We compute and store the ionization data due to the background flux ("newfluxGCRbackground") in the first column of "Ions_m2"
      % Then, we compute and store the ionization data due to the supernova flux ("newfluxdat") in the following columns for each time value
        
        Ions_m2 = zeros(Pmax,tmax+1) ;
        Ions_cm2 = zeros(Pmax,tmax+1) ;
        Ions_cm3 = zeros(Pmax,tmax+1) ;

      % Using the Simpson's Rule to numerically approximate the integral for "Ions" [ Ions_m2 = ∫(flux * ions)dE ]
        t = 1 ;
        
          for j = 1:1:Pmax

            dE = 0 ;
            I0 = 0 ;
            I1 = 0 ;
            I2 = 0 ;

            for i = 2:2:IonsEmax % (Recall, IonsEmax = fluxEmax = InterpBinNum)

              jj = Pmax-(j-1) ; % (used to flip the ionization data from sky-down to ground-up)

              dE = (Ebins(i+1)-Ebins(i-1))/2 ;
                % "dE" is multiplied within loop to compensate for not having equal intervals/steps between data points
              I0 = I0 + newfluxGCRbackground(i-1)  *newDIons(i-1,jj) *dE ; %(m2*sr*s*GeV)^(-1) * %(ion pairs) * %GeV
              I1 = I1 + newfluxGCRbackground( i )  *newDIons( i ,jj) *dE ;
              I2 = I2 + newfluxGCRbackground(i+1)  *newDIons(i+1,jj) *dE ;
              
            endfor

            Ions_m2(j,1) = (1/3)*(I0+(4*I1)+I2) ; % Simpson's Rule

            % Converting to proper units
              Ions_cm2(j,1) = Ions_m2(j,1) * (2*pi) * (10^-4) ; % (ionPairs)/(cm2*s)
                % (2pi for steradians), (10^-4 for "per m2" to "per cm2")
              Ions_cm3(j,1) = (Ions_cm2(j,1)/atmoDepth(j)) * atmoDensity(j) ; %(ionPairs)/(cm3*s)  
                %(converting to per cm3)

          endfor

        for t = 2:1:tmax+1

          for j = 1:1:Pmax

            dE = 0 ;
            I0 = 0 ;
            I1 = 0 ;
            I2 = 0 ;

            for i = 2:2:fluxEmax % (Recall, IonsEmax = fluxEmax = InterpBinNum)

              jj = Pmax-(j-1) ; % (used to flip the ionization data from sky-down to ground-up)

              dE = (Ebins(i+1)-Ebins(i-1))/2 ;
                % "dE" is multiplied within loop to compensate for not having equal intervals/steps between data points
              I0 = I0 + newfluxdat(i-1,t-1) *newDIons(i-1,jj) *dE ; %(ionPairs)/(m2*s*sr)
              I1 = I1 + newfluxdat( i ,t-1) *newDIons( i ,jj) *dE ;
              I2 = I2 + newfluxdat(i+1,t-1) *newDIons(i+1,jj) *dE ;

            endfor

            Ions_m2(j,t) = (1/3)*(I0+(4*I1)+I2) ; % Simpson's Rule

            % Converting to proper units
              Ions_cm2(j,t) = Ions_m2(j,t) * (2*pi) * (10^-4) ; % (ionPairs)/(cm2*s)
                % (2pi for steradians), (10^-4 for "per m2" to "per cm2")
              Ions_cm3(j,t) = (Ions_cm2(j,t)/atmoDepth(j)) * atmoDensity(j) ; %(ionPairs)/(cm3*s)  
                %(converting to per cm3)
          
          endfor
        endfor
        %Ions

   % Ionization Plotting
    
      if((ionsPlot) == true)
      disp("  Ionizaition Plot: ")
                      
        figure(8)
        
        hold on
        for i = 1:1:tmax+1
          semilogx(Ions_cm3(:,i), atmoAlt, "color", colors(i,:),"linewidth", 1)
        endfor
        hold off
        
        title(IonsPlotTitle)
        axis([10^(-5),10^(5),0,60])
          xlabel("Ionization [(ion pairs)/(cm^3*s)]")
          ylabel("Altitude [km]")
            set(gca, "ytick", 0:10:60)
            set(gca,'XMinorTick','on','YMinorTick','on')
          grid on
        legend(labels, "location", "northeastoutside")
  
        cd("Plots")  
        figfile = "Plot - SNCR_Ions" ;
        print(figfile,"-dpng") ;
        cd ..


      endif %plotting

      

  % %Read-in 23-Feb-1956 Ionization Rate Data
  %   if Ions1956flag = true 
  %     ions1956 = load('1956ions-cm3-s.txt')
  %   endif
  %   %disp(ions1956
  %
  %
  %
  % %%%plot the altitude of maximum ionzation vs. primary energy:
  %   if(plotMaxIonAlt eq 'yes') then begin
  %     % ..........
  %   endif %plotMaxIonAlt


  % ----------------------- Saving Ionization Data ------------------------

    if (ionsOutput == true)
    disp("Ionization Output: ")

      cd("OutputFiles")
    
      ionsFileOutput = zeros(Pmax+1,(tmax+1)+1) ;
      
      for i = 1:1:Pmax
        for j = 1:1:tmax
          ionsFileOutput(i+1,1) = atmoAlt(i) ;
         %ionsFileOutput(1,j+1) = "GCR background"
          ionsFileOutput(1,j+2) = tdat_yrs(j) ;
        endfor
        for j = 1:1:tmax+1
          ionsFileOutput(i+1,j+1) = Ions_cm3(i,j) ;
        endfor
      endfor
      
      dlmwrite(sncrionsCSV, ionsFileOutput, ",") ;
      save ("-ascii", sncrionsDAT, "Ions_cm3") ;
      
     %for j = 1 (for the GCR data)
        ionsFileOutput2 = zeros(Pmax,1) ;
        ionsFileOutput3 = zeros(Pmax,1) ;
        ionsFileOutput4 = zeros(Pmax,1) ;

        for i = 1:1:Pmax
          ionsFileOutput2(i,1) = Ions_cm2(i,1) ;
          ionsFileOutput3(i,1) = Ions_cm3(i,1) ;
          save ("-ascii", sncrionsDAT_GCR_cm2, "ionsFileOutput2") ;
          save ("-ascii", sncrionsDAT_GCR_cm3, "ionsFileOutput3") ;
        endfor
      
      for j = 2:1:tmax+1 %(start at 2 to avoid the GCR data)
        ionsFileOutput2 = zeros(Pmax,1) ;
        ionsFileOutput3 = zeros(Pmax,1) ;

        for i = 1:1:Pmax
          ionsFileOutput2(i,1) = Ions_cm2(i,j) ;
          ionsFileOutput3(i,1) = Ions_cm3(i,j) ;
          figfile2 = sprintf(sncrionsDAT_SNCR_cm2, tdat_yrs(j-1));
          figfile3 = sprintf(sncrionsDAT_SNCR_cm3, tdat_yrs(j-1));
          save ("-ascii", figfile2, "ionsFileOutput2") ;
          save ("-ascii", figfile3, "ionsFileOutput3") ;
        endfor
      endfor
    
      cd ..
        
    endif %ionsOutput
        
  endif %ions