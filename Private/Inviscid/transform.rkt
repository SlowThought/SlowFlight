#lang racket/base
#| This file is part of zfoil, written by Patrick King. It is covered by LGPL 3+. Mention my name. I'll be happy.

   A transform between z space (where cylinders live) and w space ("warped" space, where airfoils live), is developed
   by fitting a Joukowsky (Joukowski/Zhukovsky) transform to a set of points, and then applying a correction. The derivative
   of this transform is used to map the complex potential about a circle to the potential about the airfoil represented by
   the set of input points.
|#

(require racket/math              ; for π
         racket/vector            ; for vector-map
         "./potential.rkt"        ; for circle generation
         "../Utility/powell.rkt") ; for optimization routines

;;; Fit transform to set of input points. ws are values of w in "warped", airfoil space

(define (fit-transform-to-points ws θs j0 c0)
  (define n (vector-length ws))
  (define n-1 (sub1 n))
  ;; Find trailing edge - to account for blunt trailing edge, or refined calculations using displacement boundary
  ;; layer thicknesses, we take the te to be the average of the first and last ws.
  (define wte (/ (+ (vector-ref ws 0)
                    (vector-ref ws n-1))
                 2.))
; Following estimates belong elsewhere in the call chain
;  ;; Estimate leading edge (it can be a really gross estimate, wouldn't use it for chord/coefficient calculation)
;  (define wle (vector-ref ws (quotient n 2)))
;  ;; guess inititial values for j0, c0 (the center of the Joukowsky transform, and the center of the circle to be mapped)
;  (define j0 (/ (+ wle wte) 2.))
;  (define c0 (+ j0 -0.01+0.01i)) ; somewhat reasonable for airfoils of chord 1, should work in any case
  ; compute dependent parameters
  ; create (J z), (C θ)
  (define-values (J b2 zte)(make-J-transform j0 wte))
  (define a (- zte c0))
  (define C (make-circle c0 a))
  ; Improve θ estimates
  ; may need more than one pass
  (for [(i (in-range 1 (sub1 n-1)))]
    (vector-set! θs i (find-best-θ (λ (t) (J (C t)))
                                   (vector-ref ws i)
                                   (vector-ref θs (sub1 i))
                                   (vector-ref θs i)
                                   (vector-ref θs (add1 i)))))
  ; The exact transform w = T(z) = J(z)+K(z). Find K for specific ws
  (define Ks (vector-map (λ (w θ) (- w (J (C θ)))) ws θs))
  ; Return what we've learned
  (values wte J b2 zte a C θs Ks))

; given j0, wte, compute suitable Joukowsky transform and useful parameters
(define (make-J-transform j0 wte)
  (define b2 (/ (* (- wte j0)(- wte j0)) 4))
  (define zte (/ (+ j0 wte) 2))
  (values (λ (z) (+ z (/ b2 (- z j0)))) b2 zte))

; find θ such that f(θ) closest to w, given 3 guesses at θ
(define (find-best-θ f w t1 t2 t3)
  (define g (λ (t) (norm2 (- w (f t)))))
  (find-minimum (g t1)(g t2)(g t3)
                t1 t2 t3))

; fast version of (* (magnitude z)(magnitude z))
(define (norm2 z)
  (define r (real-part z))
  (define i (imag-part z))
  (+ (* r r)(* i i)))

; test code
(define ws #( 1. 0.8+0.04i 0.6+0.08i 0.4+0.11i 0.3+0.12i 0.2+0.11i 0.1+0.09i
              0.0+0.06i 0.1+0.03i 0.2+0.01i 0.3 0.4 0.6 0.8 1.0))
(define θs (build-vector 15 (λ (i) (/(* i 2 pi) 14))))
(define j0 0.5+0.03i)
(define c0 0.49+0.02i)
 
(define-values (wte J b2 zte a C new-θs Ks)
  (fit-transform-to-points ws θs j0 c0))