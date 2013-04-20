#lang racket/base
#| "circles.rkt" defines stuff related to geometry of circles

   a is a parameter of the Joukowski transform, and the complex point through which all
   Joukowski airfoil generating circles (JAGCs) pass.

   e is for "eccentricity", an old term. It defines the center of a JAGC.

   n is the number of points on a JAGC. |#
(provide radius JAGC)
(require racket/math)

;; Find the radius of a JAGC, aka the magnitude of the vector from e to a
(define(radius a e)
  (magnitude(- a e)))

;; Generate equally spaced points on a circle that could be used to generate a Joukowski airfoil
;; The last point will equal the first point will equal a (within rounding errors)
;; For a well defined leading edge, it is best to use an odd value for n.
(define(JAGC a e n)
  (let[(r(radius a e))]
    (for/vector
        [(angle(in-range(angle(- a e))7.(/ pi 0.5(sub1 n)))); Upper limit is never reached
         (i(in-range 0 n))]
      (+ e(* r(exp(* 0+i angle)))))))
(module+ test
  (printf"circle.rkt: Tests running.~n")
  (require"test-private.rkt")
  (let*[(n 21)
        (C(JAGC 0.25 -0.01 n))]
    (check-eq? n (vector-length C))
    (check-complex-=(vector-ref C 0)
                    (vector-ref C (sub1 n))
                    epsilon)))
(module+ main
  (require(submod ".." test)))

   