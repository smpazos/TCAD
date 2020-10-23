*-- #setdep @previous@

*-- @tmodel@

#define _SweepSteps_     "@Sweep_Steps@"
#define _INTENSITY_      @Intensity@
#define _WAVELENGTH_      @Lambda@

#if "@tmodel@" == "DD"
#define _TRANSMOD_     * DriftDiffusion
#define _DFORCE_       GradQuasiFermi
#define _EQUATIONSET_  Poisson Electron Hole
#elif "@tmodel@" == "HD"
#define _TRANSMOD_     Hydrodynamic
#define _DFORCE_       CarrierTempDrive
#define _EQUATIONSET_  Poisson Electron Hole eTemperature hTemperature Temperature
#elif "@tmodel@" == "Thermo"
#define _TRANSMOD_     Thermodynamic 
#define _DFORCE_       GradQuasiFermi
#define _EQUATIONSET_  Poisson Electron Hole Temperature
#endif

#-- quantum correction
#if "@QC@" == "DG"
#define _QC_ eQuantumPotential
#else
#define _QC_
#endif

#define _VREVERSE_ @Vr@
#define _VNWELL_ @VNwell@

*------------------Start of Sections--------------------*

File{
   Grid      = "@tdr@"
   *IlluminationSpectrum = "@pwd@/am15g.txt"
   Plot      = "@tdrdat@"
   *Parameter = "@parameter@"
   Current   = "@plot@"
   Output    = "@log@"
   *SpectralPlot         = "n@node@_spec"
   OpticsOutput = "n@node@_optics"
}

Electrode{
   { Name="pd1"    Voltage=_VREVERSE_ }
   { Name="pd2"    Voltage=_VREVERSE_ }
   { Name="pd3"    Voltage=_VREVERSE_ }
   { Name="pd4"    Voltage=_VREVERSE_ }
   { Name="pd5"    Voltage=_VREVERSE_ }
   { Name="pd6"    Voltage=_VREVERSE_ }
   { Name="pd7"    Voltage=_VREVERSE_ }
   { Name="pd8"    Voltage=_VREVERSE_ }
   { Name="nwell"     Voltage=_VNWELL_ }
   { Name="substrate"     Voltage=0.0 }
}

Physics{
*--Basic Physics
   AreaFactor = 1
*   HeatCapacity(TempDep)
   _TRANSMOD_  
   _QC_
   Fermi
   EffectiveIntrinsicDensity( OldSlotboom )     

*   eQuantumPotential(AutoOrientation density)
*   hQuantumPotential(AutoOrientation density)

*--Transport Models   
   Mobility(
*      CarrierCarrierScattering
      DopingDep
*      HighFieldSaturation
*      Enormal( 
*         Coulomb2D
*         IALmob
*      )
   )

*--Gen-Rec
   Recombination(
      SRH( DopingDep )
*      Auger
*      Avalanche( _DFORCE_ )
   )

*--HeavyIon Impact
*   HeavyIon (
*      Direction=(1,0)
*      Location=(0,0)
*      Time=5e-12
*      Length=17
*      Wt_hi=0.02
*      LET_f=0.0414
*      Gaussian
*      PicoCoulomb
*   )

Optics (
    ComplexRefractiveIndex (WavelengthDep(Real Imag))
    OpticalGeneration (
      QuantumYield(StepFunction(EffectiveBandgap)) * generated carriers/photon, default: 1
      ComputeFromMonochromaticSource
    ) * end OpticalGeneration
    
    Excitation (
      Intensity =  _INTENSITY_          *W/cm²
      Wavelength = _WAVELENGTH_     *um
      *Theta =
      *Phi =
      *PolarizationAngle [ deg ] or Polarization
      fromTop
      *
      Window ( 
        * origin is not required, as fromTop automatically shifts the 
        * illumination window to the top of the structure.
        # Origin = (0,-0.075)  
	*Line(x1=-1.9 x2=1.9)
        * improve accuracy due to mismatch of illumination window width & grid.
        WeightedAPDIntegration(
          * to increase accuracy, for layers thicker 100um, set > 1e4 ... 1e5
          # NumberOfCellsPerLayer = 50000 
        ) 
      ) * end Window
    ) * end Excitation

    OpticalSolver (
      TMM (
        LayerStackExtraction ()
        *IntensityPattern = Envelope
      )*end TMM
    )*end OpticalSolver
  )*end Optics
   
}*endPhysics

*Thermode{ { Name="drain" Temperature=300 SurfaceResistance=5e-4 } }
*Thermode{ { Name="source" Temperature=300 SurfaceResistance=5e-4 } }
*Thermode{ { Name="substrate" Temperature=300 SurfaceResistance=5e-4 } }

Plot{
*--Density and Currents, etc
   eDensity hDensity
   TotalCurrent/Vector eCurrent/Vector hCurrent/Vector
   eMobility hMobility
   eVelocity hVelocity
*   eQuasiFermi hQuasiFermi

*--Temperature 
*   eTemperature Temperature hTemperature

*--Fields and charges
   ElectricField/Vector Potential SpaceCharge

*--Doping Profiles
   Doping *DonorConcentration AcceptorConcentration

*--Generation/Recombination
   SRH *Auger *Band2Band
   eLifeTime hLifeTime 
*   AvalancheGeneration eAvalancheGeneration hAvalancheGeneration

*--Driving forces
   eGradQuasiFermi/Vector hGradQuasiFermi/Vector
   eEparallel hEparallel eENormal hENormal
*   eIonIntegral hIonIntegral MeanIonIntegral

*--Band structure/Composition
*   BandGap 
*   BandGapNarrowing
*   Affinity
*   ConductionBand ValenceBand
*   eQuantumPotential

*--NonLocal Meshes
*   NonLocal

*--Heavy Ion Info
*   HeavyIonChargeDensity
*   HeavyIonGeneration

*--Optics
   OpticalIntensity
   OpticalGeneration
   ComplexRefractiveIndex
}*endPlot

*CurrentPlot{
*   AvalancheGeneration( Integrate(Semiconductor) )
*}

Math {
*--Solver
   Extrapolate
   Derivatives
*   Avalderivatives
   Iterations= 15
   Notdamped= 100
   Method= Blocked
   SubMethod= Pardiso
   Number_of_Threads= 6
   StackSize = 200000000

*--Convergence Aids  
   Digits= 7            *--Relative Error=10^-Digits
*   CDensityMin= 10      *--Low Density zones
*   ErrRef(Temperature) = 30  *--Temperature reference for relative tolerated error of Temperature eq.
   ErrRef(Poisson) = 0.00258 *--Energy level reference for relative tolerated error of poisson's eq.
   ErrRef(Electron)=1e8  *--Low Density zones
   ErrRef(Hole)=1e8      *--Low Density zones
   Error(electron)=1e8  *--Low Density zones
   Error(hole)=1e8      *--Low Density zones
*   eDrForceRefDens=1e10   *--Low Density zones: high field saturation convergence
*   hDrForceRefDens=1e10   *--Low Density zones: high field saturation convergence
*   RelTermMinDensity= 1e4  *--Unrealistically High Carrier Temperature
*   RelTermMinDensityZero= 1e8  *--Unrealistically High Carrier Temperature
   CheckRhsAfterUpdate 
  
   Wallclock *--Show times in wall clock units instead of cpu time
*   Transient=BE

*--Improved Integration (for HeavyIon)
*   RecBoxIntegr
*--Ionization Integrals (for Approximate Breakdown Analysis)
*   ComputeIonizationIntegrals
*   BreakAtIonIntegral
*   AvalPostProcessing

*--Force Pseudo 3d simulation
*   _CYLINDRICAL_
}*endMath


Solve {
   *- Build-up of initial solution:
   NewCurrentPrefix="init_"
   Optics
   Coupled(Iterations=100){ Poisson _QC_ }
*   Coupled{ _EQUATIONSET_ }
   Coupled{ Poisson Electron }
   NewCurrentPrefix="results_"
   Coupled{ Poisson Electron Hole}

   *- Bias drain and gate to target bias
* NewCurrentPrefix="pd1_"   
*Quasistationary(
*      InitialStep=0.01 MinStep=1e-7 MaxStep=0.1 Increment=2
*      Goal{ Name="pd1"  Voltage= _VREVERSE_	  }
*      *Plot {Range = (0 1) Intervals=1}
*   ){ Coupled {  _EQUATIONSET_ } 
*         CurrentPlot(Time=(Range=(0 1) Intervals=1)) 
*    }

*NewCurrentPrefix="pd2_"
*Quasistationary(
*      InitialStep=0.01 MinStep=1e-7 MaxStep=0.1 Increment=2
*      Goal{ Name="pd2" Voltage= _VREVERSE_  }
*      *Plot {Range = (0 1) Intervals=1}
*   ){ Coupled {  _EQUATIONSET_ } 
*         CurrentPlot(Time=(Range=(0 1) Intervals=1)) 
*    }
*   Save(FilePrefix="InitialState_")

*NewCurrentPrefix="pd3_"
*Quasistationary(
*      InitialStep=0.01 MinStep=1e-7 MaxStep=0.1 Increment=2
*      Goal{ Name="pd3" Voltage= _VREVERSE_  }
*      *Plot {Range = (0 1) Intervals=1}
*   ){ Coupled {  _EQUATIONSET_ } 
*         CurrentPlot(Time=(Range=(0 1) Intervals=1)) 
*    }
*   Save(FilePrefix="InitialState_")

*   NewCurrentPrefix="IdTran1_"
*   Transient( 
*     InitialTime= 0 FinalTime=10e-12
*      InitialStep=1e-14 MaxStep=1e-12 MinStep=1e-16
*      Plot { Range = (0 4e-12) Intervals=4 }
*      Plot { Range = (4e-12 6e-12) Intervals=10 }
*   ) { Coupled { _EQUATIONSET_ } }
     

*   NewCurrentPrefix="IdTran2_"
*   Transient( 
*      InitialTime=10e-12  FinalTime=1e-7
*      InitialStep=1e-11 MaxStep=2e-10 MinStep=1e-12
*      Plot { Range = (10e-12 1e-7) Intervals=10 } 
*   ) { Coupled { _EQUATIONSET_ } }

*   NewCurrentPrefix="IdTran3_"
*   Transient( 
*      InitialTime=1e-7  FinalTime=1e-6
*      InitialStep=1e-9 MaxStep=1e-7 MinStep=1e-11
*      Plot { Range = (1e-7 1e-5) Intervals=10 } 
*   ) { Coupled { _EQUATIONSET_ } }
}*endSolve