SUPERNOVA COSMIC RAY ATMOSPHERIC FLUX & IONIZATION readme (10/09/2020)

For background info, see publications:
  Thomas et. al, ApJ Letters, 826, 1, 2016
    https://arxiv.org/abs/1605.04926
  Melott et al., ApJ, 840, 105, 2017
    https://arxiv.org/abs/1702.04365

Created by Alex Yelland and Brian Thomas, Washburn University, 2020
  Shared under Creative Commons license:
    Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
    https://creativecommons.org/licenses/by-nc-sa/4.0/

Questions should be directed to Brian Thomas, brian.thomas@washburn.edu


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%--- Purpose of the Program ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  This code uses data from the two papers referenced above to compute the background 
  galactic cosmic ray (GCR) flux, the minimum (Case B) supernova cosmic ray (SNCR) flux, 
  and the resulting atmospheric ionization for various time periods.
    
    {Case B = diffusive case for CR transport, disordered magnetic field 
              between the cosmic ray source (the supernova) and Earth}
  
  By altering the distance between the CR source and Earth, we can examine a variety of 
  cases of the resulting flux and ionization data.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%--- Structure of the Program ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  NOTE: Before compiling, ensure that the source code "SNCR_Flux_and_Ionization.m" is in 
        the same directory as the "OutputFiles", "ReadinFiles", and "Plots" folders. Do 
        not change the name of the folders unless you change the references within the 
        source code as well.

  Main Editing Section: ------------------------------------------------------------------

    - Throughout the code, all file names are referenced using the variables in this 
    section; therefore, if wanting to change the name of any input/output file, change it 
    here.
    - The code is broken up into multiple sections that can be "turned-off" or "turned-on" 
    during compiling. By changing the values of the flags (1=on, 0=off), different 
    sections of the code can be skipped. This includes the plotting codes at the end of 
    each section.

  Cosmic Ray Flux: -----------------------------------------------------------------------
    
    Background GCR Flux:
      - Read-in the time periods and energy bin levels we want to analyze our flux and 
      ionization data at. 
      - Use energy bin date in the calculation of the background GCR flux (Haungs: Rebel 
      and Roth: 2003: Reports on Progress in Physics: vol. 66: p. 1145-1206). This 
      calculation consist of two power laws (one above the "knee" and below the "knee") 
      before it is then multiplied by E^2 to scale the data for the comparison to the SNCR
      flux.

    Supernova CR Flux:
      - Using data from Melott et al., ApJ, 840, 105, 2017, the supernova's cosmic ray 
      flux is calculated at the same time periods and energy bin values as the 
      background GCR flux.
    
    Flux Output:
      - The SNCR flux is saved as ".csv" with the first row labeling the time period
      and the first column labeling the energy bin value. 
      - The GCR flux is saved into two ".dat" files (one with only GCR flux data, and the 
      other with the GCR flux data and its associated energy bin values)

  Ions: ----------------------------------------------------------------------------------
    
    Read-in CR Flux Data:
      - This section of the code is only compiled when the Cosmic Ray Flux calculation section
      is skipped (flux=0). In order to calculate the atmospheric ionization, we need some 
      form of CR flux. So, this section of the code allows you to read-in GCR flux values
      and SNCR flux values rather than computing within the code.
    
    Read-in Dimitra Atri's Ionization Data:
      - Dimitra's (symbolized with a "D") ionization data represents the ionization of a 
      single cosmic ray (a proton) entering the atmosphere causing the ionizing "shower"
      effect.  Calculated from CORSIKA runs, see Atri, D. et al., 2010, JCAP, 2010, 008
      - Read-in the single CR ionization data along with its specified energy bin 
      values and altitude levels.
    
    Interpolation:
      - Choose a larger, odd integer for the interpolation parameter (InterpBinNum). 
      - At this value, interpolate to a new higher resolution energy bin array (Ebins) 
      with upper and lower bounds determined by the SNCR flux energy bins (Edat) and 
      Dimitra's energy bins (DEdat). 
      - Repeat this process for the SNCR flux data and Dimitra's single CR ionization
      data using a logarithmic interpolation (see "log_interpol.pdf"). 
      - The higher resolution energy bin values (Ebins) is used to recalculate the GCR
      flux. 
      - After interpolation, all read-in data and calculated data on the same 
      resolution scale.
    
    Atmospheric Ionization Calculation:
      - Use the Simpson's Rule to numerically approximate the overall atmospheric 
      ionization for the background galactic cosmic rays and supernova cosmic rays at 
      different altitudes and time periods. 
      - The integral for this calculation is "Ions = ∫(flux * Dions)*dE"
      - Once calculated, use the atmospheric pressure data to calculate the ionization 
      volume density.
    
    Ionization Output:
      - The atmospheric ionization data is saved as ".csv" with the first row labeling the
      time period and the first column labeling altitude. The GCR ionization is saved in 
      the second column. 
      - The ionization data is then broken apart into the GCR ionization and the 
      individual time periods' ionizations and saved as ".dat" files.
