; UCS coordinates
(sde:set-process-up-direction 1)

;----------------------------------------------------------------------
; Setting parameters
; - lateral
; 18 u de ancho de difusión y 1 um de separación entre bordes de difusión.

(define Ddiff @Diff_Spacing@)           	;diffusion minimum drawn spacing ON C5
(define Hepi @Epi_Depth@)			;[um] EPI depth
(define Diff_count 	8)			;Number of diffusions in the structure
(define Rdiff 8.25)				;Diffusion half width
(define Totdiff	(* 2 Rdiff))			;Diffusion total width
(define Nwellwidth 18.5)                        ; [um] Nwell width
(define Ltot (+ (+ (* (+ Totdiff Ddiff) Diff_count) Nwellwidth) Rdiff)) ; [um] Lateral extend total
(define HLocos 0.29)
(define fillet-radius (/ HLocos 2))          ; [um] Rounding radius

; - layers
(define Hsub 40)                        ; [um] Substrate thickness
(define Nwelldepth 3)                        ; [um] N-well depth
(define OxyNithk 6)                        ; [um] Substrate thickness
(define ILox 5e-4)                     ; [um] Sacrificial contact defining Layer thickness

; - pn junction resolution
(define Gpn 0.005)                  ; [um]


; - Substrate doping level
;(define Pepi 	2e16)                 ; [1/cm3]
(define Pepi 	@Epi_Doping@)                 ; [1/cm3]
;(define Pepi 	1.5e17)                 ; [1/cm3]
(define Psub 	1e19)                 ; [1/cm3]
(define Nplus  	5e18)                 ; [1/cm3]
;(define Nplus  	1e19)                 ; [1/cm3]
(define Xj  0.25)                    ; [um] Junction depth
(define Pchstop 	1e19)                 ; [1/cm3]
(define Nwell  	1e17)                 ; [1/cm3]

; dopants
(define SD-dopant        "ArsenicActiveConcentration")
(define substrate-dopant "BoronActiveConcentration")

;----------------------------------------------------------------------
; Derived quantities
(define Ymax Ltot)
;(define Ydiff (/ Ymax 2))	
(define Ydiff 0)	

;(define Ydmax (+ Ydiff Rdiff))	;y-Coordinate for diffusion limit right
(define Ydmin (+ Ddiff Totdiff))	;y-Coordinate for diffusion limit left
(define Ydmax (+ Ydmin Totdiff))	;y-Coordinate for diffusion limit right

(define Xsub Hsub)		;x-Coordinate for substrate bottom limit
(define Xepi Hepi)		;x-Coordinate for EPI bottom limit
(define XILox (- ILox))		;x-Coordinate for scarficial contact defining layer 
(define XOxyNi (- OxyNithk))	;x-Coordinate for scarficial contact defining layer 
(define XLocos (/ HLocos 2))	;x-Coordinate for bottom of LOCOS structure 
(define XChstop (- Xj XLocos))  ;Depth of the Channel stop implants

(define COUNT (list 1 2 3 4 5 6 7))
(define COUNTS (list "2" "3" "4" "5" "6" "7" "8"))


;----------------------------------------------------------------------
; procedure definitions
;

; convenience procedure for specifying 2-d coordinates
; this function defined to avoid need to specify 0 coordinate for z
;
(define (position-2d x y)
  (position x y 0)
)

;----------------------------------------------------------------------
; Create the physical structure.  This is the half-FET.  The full
; structure is symmetric and will be generated by reflecting the
; half structure at the end of the meshing process
;

; Creating epi-substrate
(define substrate-region
  (sdegeo:create-rectangle
    (position-2d Xepi (* 0 Ymax))
    (position-2d 0 Ymax)
    "Silicon"
    "R.Epi"
  )
)

; Creating substrate
(define substrate-region
  (sdegeo:create-rectangle
    (position-2d Xsub (* 0 Ymax))
    (position-2d Xepi Ymax)
    "Silicon"
    "R.Substrate"
  )
)

; Creating FOX nitride
(define IL-oxide
  (sdegeo:create-rectangle
    (position-2d 0 0)
    (position-2d XOxyNi Ymax)        ;-6 um oxynitride (SEM en Uddin2013)
    "Nitride"
    "R.FOX"
  )
)


; Creating simplified LOCOS oxide
(define IL-oxide
  (sdegeo:create-rectangle
    (position-2d 0 (- 0 5))
    (position-2d XLocos (- Ddiff 5))        
    "SiO2"
    "R.LOCOS"
  )
)
(sdegeo:fillet-2d
  (find-vertex-id (position-2d XLocos (- 0 5)))
  fillet-radius
)
(sdegeo:fillet-2d
  (find-vertex-id (position-2d XLocos (- Ddiff 5)))
  fillet-radius
)
(sdegeo:move-2d-regions (find-body-id (position-2d (/ XLocos 2) (- 0 4.9))) (gvector 0 (+ 5 Totdiff) 0))

(for-each
	(lambda (COUNT)
		(begin
			; Creating simplified LOCOS oxide
			(define IL-oxide
			  (sdegeo:create-rectangle
			    (position-2d 0 (- 0 5))
			    (position-2d XLocos (- Ddiff 5))        
			    "SiO2"
			    "R.LOCOS"
			  )
			)
			(sdegeo:fillet-2d
			  (find-vertex-id (position-2d XLocos (- 0 5)))
			  fillet-radius
			)
			(sdegeo:fillet-2d
			  (find-vertex-id (position-2d XLocos (- Ddiff 5)))
			  fillet-radius
			)
			(sdegeo:move-2d-regions (find-body-id (position-2d (/ XLocos 2) (- 0 4.9))) (gvector 0 (+ 5 (+ Totdiff (* Ydmin COUNT))) 0))
		)
	)COUNT
)

(define IL-oxide
  (sdegeo:create-rectangle
    (position-2d 0 (- 0 5))
    (position-2d XLocos (- 1.5 5))        
    "SiO2"
    "R.LOCOS"
  )
)
(sdegeo:fillet-2d
  (find-vertex-id (position-2d XLocos (- 0 5)))
  fillet-radius
)
(sdegeo:fillet-2d
  (find-vertex-id (position-2d XLocos (- 1.5 5)))
  fillet-radius
)
(sdegeo:move-2d-regions (find-body-id (position-2d (/ XLocos 2) (- 0 4.9))) (gvector 0 (+ 5 (+ Totdiff (* Ydmin (- Diff_count 1)))) 0))

;**********Sacrificial shapes for contact definition************
;***************************************************************
; Creating Interfacial Layer oxide
(define IL-oxide
  (sdegeo:create-rectangle
    (position-2d 0 (- Totdiff 4))
    (position-2d XILox (- Totdiff 2))
    "SiO2"
    "R.ILox"
  )
)
;Delete the axiliary oxide layer to define contact region
(sdegeo:delete-region (list (find-region-id "R.ILox")))

(for-each
	(lambda (COUNT)
		(begin
			; Creating Interfacial Layer oxide
			(define IL-oxide
			  (sdegeo:create-rectangle
			    (position-2d 0 (- (+ Totdiff (* Ydmin COUNT)) 4))
			    (position-2d XILox (- (+ Totdiff (* Ydmin COUNT)) 2))
			    "SiO2"
			    "R.ILox"
			  )
			)
			;Delete the axiliary oxide layer to define contact region
			(sdegeo:delete-region (list (find-region-id "R.ILox")))
		)
	)COUNT
)

(define IL-oxide
  (sdegeo:create-rectangle
    (position-2d 0 (- (- Ltot Rdiff) 2))
    (position-2d XILox (- (- Ltot Rdiff) 4))
    "SiO2"
    "R.ILox"
  )
)
;Delete the axiliary oxide layer to define contact region
(sdegeo:delete-region (list (find-region-id "R.ILox")))

(define IL-oxide
  (sdegeo:create-rectangle
    (position-2d 0 (- Ltot 0.5))
    (position-2d XILox (- Ltot 1.5))
    "SiO2"
    "R.ILox"
  )
)
;Delete the axiliary oxide layer to define contact region
(sdegeo:delete-region (list (find-region-id "R.ILox")))
;----------------------------------------------------------------------
; Contact formation: define then assign
;
(sdegeo:define-contact-set
  "pd1"
)

(sdegeo:define-contact-set
  "substrate"
)

(sdegeo:define-contact-set
  "nwell"
)

(sdegeo:define-2d-contact
 (find-edge-id (position-2d 0 (- Totdiff 3)))
 "pd1"
)

(for-each
	(lambda (COUNT COUNTS)
		(begin
			(define CONTNAME (string-append "pd" COUNTS))
			(sdegeo:define-contact-set CONTNAME)
			(sdegeo:define-2d-contact
			 (find-edge-id (position-2d 0 (- (+ Totdiff (* Ydmin COUNT)) 3)))
			 CONTNAME
			)
		)
	)COUNT COUNTS
)

(sdegeo:define-2d-contact
 (find-edge-id (position-2d 0 (- (- Ltot Rdiff) 3)))
 "nwell"
)

(sdegeo:define-2d-contact
 ;(find-edge-id (position-2d Xsub 5e-4))
 (find-edge-id (position-2d 0 (- Ltot 1)))
 "substrate"
)


;----------------------------------------------------------------------
; Doping profiles:
;

;;; - Epi (Constant Profile)
(sdedr:define-constant-profile
  "Const.Epi"
  substrate-dopant
  Pepi
)

(sdedr:define-constant-profile-region
  "PlaceCD.Epi"
  "Const.Epi"
  "R.Epi"
)

;;; - Substrate (Constant Profile)
(sdedr:define-constant-profile
  "Const.Substrate"
  substrate-dopant
  Psub
)

(sdedr:define-constant-profile-region
  "PlaceCD.Substrate"
  "Const.Substrate"
  "R.Substrate"
)

; n+ Photodiode implant: definition of profile
(sdedr:define-gaussian-profile
  "Impl.SDprof"
  SD-dopant
  "PeakPos" 0 
  "PeakVal" Nplus
  "ValueAtDepth" Pepi
  "Depth" Xj
  "Gauss"
  "Factor" 0.2
)

; p+ Channel Stop implant: definition of profile
(sdedr:define-gaussian-profile
  "Impl.ChStop"
  substrate-dopant
  "PeakPos" 0 
  "PeakVal" Pchstop
  "ValueAtDepth" 1e18
  "Depth" (/ Xj 1.5)
  "Gauss"
  "Factor" 0.002
)

; n-well ring implant: definition of profile
(sdedr:define-gaussian-profile
  "Impl.Nwell"
  SD-dopant
  "PeakPos" 0 
  "PeakVal" Nwell
  "ValueAtDepth" Pepi
  "Depth" Nwelldepth
  "Gauss"
  "Factor" 0.2
)

; Source/Drain base line definition: direction matters
(sdedr:define-refinement-window "BaseLine.Pd1" "Line"
 (position-2d 0 Totdiff)
 (position-2d 0 0)
)

; Source/Drain implant: placement of profile to baseline
(sdedr:define-analytical-profile-placement
  "Impl.Pd1"
  "Impl.SDprof"
  "BaseLine.Pd1"
  "Positive"
  "NoReplace"
  "Eval"
)

; Channel stop base line definition: direction matters
(sdedr:define-refinement-window "BaseLine.ChStopL" "Line"
 (position-2d 0 Ydmin)
 (position-2d 0 Totdiff)
)

; Source/Drain implant: placement of profile to baseline
(sdedr:define-analytical-profile-placement
  "Impl.ChStopL"
  "Impl.ChStop"
  "BaseLine.ChStopL"
  "Positive"
  "NoReplace"
  "Eval"
)

; Nwell base line definition: direction matters
(sdedr:define-refinement-window "BaseLine.Nwell" "Line"
 (position-2d 0 (+ (- Ltot Rdiff) 0.5))
 (position-2d 0 (- (+ (- Ltot Rdiff) 0.5) Nwellwidth))
)

; Nwell implant: placement of profile to baseline
(sdedr:define-analytical-profile-placement
  "Impl.Nwellplace"
  "Impl.Nwell"
  "BaseLine.Nwell"
  "Positive"
  "NoReplace"
  "Eval"
)

(for-each
	(lambda (COUNT COUNTS)
		(begin
			(define IMPNAME (string-append "Impl.Pd" COUNTS))
			(define CSIMPNAME (string-append "Impl.ChStop" COUNTS))
			(define BLINENAME (string-append "BaseLine.Pd" COUNTS))
			(define CSBLINENAME (string-append "BaseLine.ChStop" COUNTS))
			;********************
			;Photodiode implants
			;********************
			; Source/Drain base line definition: direction matters
			(sdedr:define-refinement-window BLINENAME "Line"
			 (position-2d 0 (+ Totdiff (* Ydmin COUNT)))
			 (position-2d 0 (* Ydmin COUNT))
			)

			; Source/Drain implant: placement of profile to baseline
			(sdedr:define-analytical-profile-placement
			  IMPNAME
			  "Impl.SDprof"
			  BLINENAME
			  "Positive"
			  "NoReplace"
			  "Eval"
			)

			;********************
			;Chanel stop implants
			;********************
			; Channel stop base line definition: direction matters
			(sdedr:define-refinement-window CSBLINENAME "Line"
			 (position-2d 0 (+ (+ Totdiff (* Ydmin COUNT)) Ddiff))			 
			 (position-2d 0 (+ Totdiff (* Ydmin COUNT)))
			)

			; Channel stop implant: placement of profile to baseline
			(sdedr:define-analytical-profile-placement
			  CSIMPNAME
			  "Impl.ChStop"
			  CSBLINENAME
			  "Positive"
			  "NoReplace"
			  "Eval"
			)
		)
	)COUNT COUNTS
)

(define CSIMPNAME (string-append "Impl.ChStop" "8"))
(define CSBLINENAME (string-append "BaseLine.ChStop" "8"))
; Channel stop base line definition: direction matters
(sdedr:define-refinement-window CSBLINENAME "Line"
 (position-2d 0 (+ (+ Totdiff (* Ydmin 7)) 1.5))			 
 (position-2d 0 (+ Totdiff (* Ydmin 7)))
)

; Channel stop implant: placement of profile to baseline
(sdedr:define-analytical-profile-placement
  CSIMPNAME
  "Impl.ChStop"
  CSBLINENAME
  "Positive"
  "NoReplace"
  "Eval"
)
;----------------------------------------------------------------------
; Meshing Specifications
;

;
; Substrate
;

; default mesh spacing: note minimum will limit junction refinement
(sdedr:define-refinement-size
   "Ref.Substrate"
   (/ Totdiff 5) (/ Hsub 16)
   (* Gpn 1) (* Gpn 1)
)

; doping refinement along junction
(sdedr:define-refinement-function
  "Ref.Substrate"
  "DopingConcentration" "MaxTransDiff" 1
)

; place the substrate refinement
(sdedr:define-refinement-region
  "RefPlace.Substrate"
  "Ref.Substrate"
  "R.Substrate"
)

; place the substrate refinement
(sdedr:define-refinement-region
  "RefPlace.Epi"
  "Ref.Substrate"
  "R.Epi"
)

;
; Active (region near the device)
;
(sdedr:define-refinement-window
  "RWin.Act"
  "Rectangle"
  (position-2d 0 0)
  (position-2d (* Xj 3) Ymax)
  ;(position-2d (+ Nwelldepth (* Xj 3)) Ymax)
)

(sdedr:define-multibox-size "MultiboxDefinition_1" (/ Xj 6) (/ Rdiff 20) (/ Xj 12) (/ Rdiff 30) 1 1)

(sdedr:define-multibox-placement "MultiboxPlacement_1" "MultiboxDefinition_1" "RWin.Act")

;----------------------------------------------------------------------
; Build Mesh
(sde:build-mesh "snmesh" "" "n@node@_half_msh")

;----------------------------------------------------------------------
; Reflect the device
; Note this makes a system call to the Sentaurus Datex Explorer
; there are other approaches, but this one guarantees the mesh
; is perfectly symmetric rather than reflecting the boundary and then
; meshing that.

;(system:command "tdx -mtt -y -M 0 -S 0 -ren pd1=pd1 pd2=pd2 n@node@_half_msh  n@node@_msh")
(system:command "tdx -mtt -y -M 0 -S 0 -same_region n@node@_half_msh  n@node@_msh")
(sde:save-model "/home/Simulations/local/sentaurus_project/2D_photoarray/sde_model")