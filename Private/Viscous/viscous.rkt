#lang racket/base
;; Private/viscous.rkt provides boundary layer, drag calculations
(provide CD)

(require "geometry.rkt"
         "inviscid.rkt"
         "quad.rkt")

(module+ test
  (require "test-private.rkt"))

(define (CD geometry alpha re #:optional Vs)
  ; Extract info from geometry structure
  (define Ss (get-S geometry))
  (define chord (geometry-user-chord geometry))
  ; Compute kinematic viscosity
  (define nu (/ re chord))
  ; If neccessary, calculate velocities
  (cond [(not Vs)
         (set! Vs (V alpha geometry))])
  ; Find stagnation point
  (define i-stag
    (let loop [(j 0)]
      (if (positive? (vector-ref Vs j))
          j
          (loop (add1 j)))))
  ; Calc upper BL backwards from stagnation point
  (define-values (upper-momentum upper-displacement upper-seperation) 
    (BL Ss Vs (sub1 i-stag) -1 nu))
  ; Calc lower BL forward from stagnation point
  (define-values (lower-momentum lower-displacement lower-seperation) 
    (BL Ss Vs i-stag +1 nu))
  ; Calc pressure drag from displacement of streamlines
  ; Combine momentum thicknesses, pressure drag to form total drag
  ; return total drag as coefficient (referenced to user chord)
  void)

(define (BL Ss Vs i0 inc nu)
  (define laminar #t)
  ; from i0 to end (either 0 or n) do
  ;    if laminar apply Thwaites Method
  ;               if transition (Michel's criteria)
  ;               or separation (Thwaite's parameter)
  ;               set laminar #f
  ;               next step
  ;    else apply White's method
  ;         if seperation return position, theta, delta*
  ;         next step
  void)

#| Boundary layer calculations are done in n,s coordinate system, where n is
   normal to airfoil surface, and s is parallel to it. The following converts
   an array of complex coordinates (as in a description of an airfoil) to an
   array of distances from the TE, with the first point defined as zero. |#
(define (Z->S Zs)
  (let [(s 0.)
        (z_old (vector-ref Zs 0))]
    (for/vector [(z Zs)]
      (set! s (+ s (magnitude (- z z_old))))
      (set! z_old z)
      s)))

(define (int-U^5 Vs Ss)
  ; Find stagnation point
  (define i-stag (let loop [(i 0)]
                   (if (< 0. (vector-ref Vs i))
                       i
                       (loop (add1 i)))))
  (interpolate (vector-ref Vs i-stag)
               (vector-ref Vs (sub1 i-stag))
               (vector-ref Ss i-stag)
               (vector-ref Ss (sub1 i-stag))
               0.))