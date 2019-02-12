#lang racket/base
#| This file is part of zfoil, written by Patrick King. It is covered by LGPL 3+. Mention my name. I'll be happy.
   Functions for describing the inviscid flow around a circle are developed using complex potentials.
|#

(require racket/math) ; for complex functions
(provide make-circle           ; complex center, complex radius -> ((function real angle) -> complex on circle)
         cylinder-parameters   ; complex radius -> magnitude radius, norm^2 radius, α0L
         Φ_cylinder            ; complex center, magnitude radius, norm^2 radius, α0L -> complex potential
         dΦ/dz_cylinder)       ; ... -> complex derivative


;; First, we define the cylinder about which we wish to determine the flow
(define (make-circle c0 a)
         (λ (θ)(+ (* a (exp (* 0+i θ))) c0)))

#| Flow about a cylinder with circulation  centered on c0. The cylinder is defined by z = ae^iθ + c0, a complex. 
   We expect a and zte to be consistent with the definitions in "joukowski.rkt" (zte = a + c0), and that these values
   are already available. |#

(define (cylinder-parameters a)
  (let*[(na (magnitude a))
        {na^2 (* na na)}
        (α0L (angle a))]
    (values na na^2 α0L)))

(define (Φ_cylinder c0 na na^2 α α0L)
  (λ(z)(+ ((Φ_Uniform_Flow (V∞ α) c0) z)
          ((Φ_Doublet α na^2 c0) z)
          ((Φ_Vortex (* 2 na (sin (- α α0L))) c0) z))))


(define (dΦ/dz_cylinder c0 na na^2 α α0L)
  (λ(z)(+ ((dΦ/dz_UF (V∞ α)) z)
          ((dΦ/dz_Dub α na^2 c0) z)
          ((dΦ/dz_Vortex (* 2 na (sin (- α α0L))) c0) z))))

(module+ test
  (printf "potential.rkt tests running.~n")
  (require plot
           "../Utility/plots.rkt")
  ; For sake of clarity in some plots, and irrelevance to the eventual application,
  ; we suppress the display of velocity inside certain regions.
  (define (suppress-f-within-r-of-c f r c) ; f complex->complex, r & c complex, returned function vector->vector
    (λ(x y)
      (let [(z (make-rectangular x y))]
        (if (< (magnitude (- z c)) r) #(0 0)
             (complex->vector (f z)))))))

(module+ test
  (define c0 1+i)
  (define a 1-i)
  (define c (make-circle c0 a))
  (plot (parametric (λ (t) (complex->vector (c t))) 0 6.29))
  (printf "~nYou should see a circle above.~n"))

;; Velocity in the far field is just uniform flow of magnitude 1
(define (V∞ α) (exp (* 0+i α))) 

;; The potential function and its derivative use the conjugate of the velocity (C-R vs continuity)
;; This is most obvious in uniform flow.

(define (Φ_Uniform_Flow Vinf c0)(λ(z)(* (conjugate Vinf) (- z c0))))
(define (dΦ/dz_UF Vinf)(λ(z)(conjugate Vinf)))

(module+ test
  (define alpha 0.3)
  (define Vinf (V∞ alpha))
  (define UF (Φ_Uniform_Flow Vinf c0))
  (define dUF (dΦ/dz_UF Vinf))
  (display
   (plot (list (parametric (λ (t) (complex->vector (c t))) 0 6.29)
               (contours (λ(x y)(let[(z (make-rectangular x y))](imag-part (UF (make-rectangular x y)))))
                         #:levels 5)
               (vector-field (λ(x y)(complex->vector (conjugate (dUF (make-rectangular x y)))))))
         #:x-min -1
         #:x-max 3
         #:y-min -1
         #:y-max 3))
  (printf "~nUniform flow, shown two ways.~n"))


;; Doublet oriented to angle of attack, of strength γ, centered on c0
(define (Φ_Doublet α γ c0)(λ(z)(/(* γ (exp (* 0+i α)))(- z c0))))
(define (dΦ/dz_Dub α γ c0)(λ(z)(/(* -1. γ (exp (* 0+i α)))(- z c0)(- z c0))))

(module+ test
  (define Dub (Φ_Doublet alpha 1 c0))
  (define dDub(dΦ/dz_Dub alpha 1 c0))
  (display (plot (list (parametric (λ (t) (complex->vector (c t))) 0 6.29)
                       (contours (λ (x y) (imag-part (Dub (make-rectangular x y))))
                                 #:levels 11)
                       (vector-field (suppress-f-within-r-of-c (λ(z)(conjugate (dDub z))) 0.5 c0)
                                     #:scale 0.1))
                 #:x-min -1
                 #:x-max 3
                 #:y-min -1
                 #:y-max 3))
  (printf "~nDoublet flow, shown two ways.~n")
  #| A uniform flow and a doublet combined should yield the flow field about a cylinder. Say that cylinder
     is of radius |a|, centered on c0, where a and c0 are complex numbers (making a complex will prove convenient later).
     γ should be |a|^2.
  |#
  (let [(na^2 (*(magnitude a)(magnitude a)))]
    (set! Dub (Φ_Doublet alpha na^2 c0))
    (set! dDub (dΦ/dz_Dub alpha na^2 c0))
    (display
     (let [(f (λ(z)(+ (UF z)(Dub z))))
           (df (λ(z)(+ (dUF z)(dDub z))))]
       (plot (list (contours (λ(x y)(imag-part (f (make-rectangular x y)))) #:levels 15)
                   (vector-field (suppress-f-within-r-of-c (λ (z)(conjugate (df z))) (magnitude a) c0)
                                 #:scale .2)
                   (parametric (λ (t) (complex->vector (c t))) 0 6.29))
           #:x-min -2
           #:x-max 4
           #:y-min -2
           #:y-max 4)))
  (printf "~nUniform flow + doublet = flow about a cylinder, no circulation.~n")))

;; Vortex centered on c0
(define (Φ_Vortex Γ c0) (λ(z) (* 0+i Γ (log (- z c0)))))
(define (dΦ/dz_Vortex Γ c0) (λ(z)(/(* 0+i Γ)(- z c0))))

(module+ test
  (define PhiV (Φ_Vortex 1 c0))
  (define dPhiV(dΦ/dz_Vortex 1 c0))
  (display (plot (list (contours (λ(x y)(imag-part (PhiV (make-rectangular x y)))))
                       (vector-field (λ(x y)(complex->vector (conjugate (dPhiV (make-rectangular x y)))))))
                 #:x-min 0
                 #:x-max 2
                 #:y-min 0
                 #:y-max 2))
  (printf "~nFlow about a vortex.~n"))



(module+ test
  (define-values (na na^2 α0L)(cylinder-parameters a))
  (printf "a is ~s~n|a| is ~s~n|a|^2 is ~s~nα0L is ~s~n" a na na^2 α0L)
  (display (plot (list (parametric (λ (t) (complex->vector (c t))) 0 6.29)
                       (contours (λ(x y)(imag-part ((Φ_cylinder c0 na na^2 alpha α0L)(make-rectangular x y))))
                                 #:levels 20)
                       (vector-field (suppress-f-within-r-of-c
                                      (λ(z)(conjugate((dΦ/dz_cylinder c0 na na^2 alpha α0L)z))) na c0)))
                 #:x-min -2
                 #:x-max 4
                 #:y-min -2
                 #:y-max 4))
  (printf "~nFlow about a cylinder with circulation.~n"))