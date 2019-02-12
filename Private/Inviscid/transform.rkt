#lang racket/base
#| This file is part of zfoil, written by Patrick King. It is covered by LGPL 3+. Mention my name. I'll be happy.

   A transform between z space (where cylinders live) and w space ("warped" space, where airfoils live), is developed
   by fitting a Joukowsky (Joukowski/Zhukovsky) transform to a set of points, and then applying a correction. The derivative
   of this transform is used to map the complex potential about a circle to the potential about the airfoil represented by
   the set of input points.
|#

(require "./potential.rkt") ; for circle generation

;;; Fit transform to set of input points. ws are values of w in "warped", airfoil space
(define (fit-transform-to-points ws)
  (define n (vector-length ws))
  (define n-1 (sub1 n))
  ;; Find trailing edge - to account for blunt trailing edge, or refined calculations using displacement boundary
  ;; layer thicknesses, we take the te to be the average of the first and last ws.
  (define wte (/ (+ (vector-ref ws 0)
                    (vector-ref ws n-1))
                 2.))
  ;; Estimate leading edge (it can be a really gross estimate, wouldn't use it for chord/coefficient calculation)
  (define wle (vector-ref ws (quotient n 2)))
  ;; guess inititial values for j0, c0 (the center of the Joukowsky transform, and the center of the circle to be mapped)
  (define j0 (/ (+ wle wte) 2))
  (define c0 (+ j0 -0.01+0.01i)) ; somewhat reasonable for airfoils of chord 1, should work in any case
  ;; given j0, c0, find an error value
   ; (let loop [(j0 j0)(c0 c0)(i 0)]
   ;    (if (= i 10)
   ;        ; return relevant parameters
   ;        (begin
              ; compute dependent parameters
              ; create (J z), (C θ)
              (define-values (J b2 zte)(make-J-transform j0 wte))
              (define a (- zte c0))
              (define C (make-circle c0 a))
     ; find points on circle that map closely to ws
     ; sum up square of distance between ws and (J(C(θs))
     ; pick new j0, c0 to reduce error
  )
  
; given j0, wte, compute suitable Joukowsky transform and useful parameters
(define (make-J-transform j0 wte)
  (define b2 (/ (* (- wte j0)(- wte j0)) 4))
  (define zte (/ (+ j0 wte) 2))
  (values (λ (z) (+ (- z j0)(/ b2 (- z j0)))) b2 zte))
