#lang racket/base
#| This file is part of zfoil, written by Patrick King. It is covered by LGPL 3+. Mention my name. I'll be happy.

   Functions for describing the inviscid flow around a circle are developed using complex potentials.
|#
(require racket/math) ; for complex functions
(module+ test
  (printf "potential.rkt tests running.~n")
  (require plot
           "../Utility/plots.rkt"))

;; Velocity in the far field is just uniform flow of magnitude 1
(define (V∞ α) (exp (* 0+i α))) 

;; The potential function and its derivative use the conjugate of the velocity (C-R vs continuity)
(define (Φ_Uniform_Flow Vinf)(λ(z)(* (conjugate Vinf) z)))
(define (dΦ/dz_UF Vinf)(λ(z)(conjugate Vinf)))

(module+ test
  (define alpha 0.3)
  (define Vinf (V∞ alpha))
  (define UF (Φ_Uniform_Flow Vinf))
  (define dUF (dΦ/dz_UF Vinf))
  (display
   (plot (list (contours (λ(x y)(imag-part (UF (make-rectangular x y)))) #:levels 5)
               (vector-field (λ(x y)(complex->vector (conjugate (dUF (make-rectangular x y)))))))
         #:x-min 0
         #:x-max 2
         #:y-min 0
         #:y-max 2))
  (printf "~nUniform flow, shown two ways.~n"))


;; Doublet oriented to angle of attack, of strength γ, centered on c0
(define (Φ_Doublet α γ c0)(λ(z)(/(* γ (exp (* 0+i α)))(- z c0))))
(define (dΦ/dz_Dub α γ c0)(λ(z)(/(* -1. γ (exp (* 0+i α)))(- z c0)(- z c0))))

(module+ test
  (define c0 1+1i)
  (define Dub (Φ_Doublet alpha 1 c0))
  (define dDub(dΦ/dz_Dub alpha 1 c0))
  (display
   (plot (list (contours (λ(x y)(imag-part (Dub (make-rectangular x y)))) #:levels 15)
               (vector-field (λ(x y)(complex->vector (conjugate (dDub (make-rectangular x y)))))
                             #:scale 0.01))
         #:x-min 0
         #:x-max 2
         #:y-min 0
         #:y-max 2))
  (printf "~nDoublet flow, shown two ways.~n")
  #| A uniform flow and a doublet combined should yield the flow field about a cylinder. Say that cylinder
     is of radius |a|, centered on c0, where a and c0 are complex number (making a complex will prove convenient later).
     γ should be |a|^2.
  |#
  (define a 0.5)
  (set! Dub (Φ_Doublet alpha (* a a) c0))
  (set! dDub (dΦ/dz_Dub alpha (* a a) c0))
  (display
   (let [(f (λ(z)(+ (UF z)(Dub z))))
         (df (λ(z)(+ (dUF z)(dDub z))))]
     (plot (list (contours (λ(x y)(imag-part (f (make-rectangular x y)))) #:levels 15)
                 (vector-field (λ(x y)
                                 ; For sake of clarity in the plot, and irrelevance to the eventual application,
                                 ; we suppress the display of velocity inside the cylinder
                                 (if (< (+ (* (sub1 x)(sub1 x))(* (sub1 y)(sub1 y)))(* a a))
                                     #(0 0)
                                     (complex->vector (conjugate (df (make-rectangular x y))))))
                               #:scale .1))
           #:x-min 0
           #:x-max 2
           #:y-min 0
           #:y-max 2)))
  (printf "~nFlow about a cylinder, no circulation.~n"))

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
