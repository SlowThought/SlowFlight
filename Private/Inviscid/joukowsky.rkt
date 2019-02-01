#lang racket/base
(provide (all-defined-out))
(require racket/math "./potential.rkt")
#| This file is part of zfoil, written by Patrick King. It is covered by LGPL 3+. Mention my name. I'll be happy.
   The Joukowsky complex transform maps circles near the origin into airfoil-like shapes, while
   leaving the far field undistorted. z space contains the circles. w (warped) space contains the
   airfoils.

   c0 is the center of the circle that is mapped to an airfoil shape (given).
   j0 is a parameter of the transform between z space and w space (given).
   wte is the location of the trailing edge of the airfoil of interest in w space (given).
   These 3 complex parameters (6 real parameters) are sufficient to determine a Joukowsky airfoil.
   We generate 3 "dependent parameters" that are often original in other constructions. |#

(define (make-dependent-parameters c0 j0 wte)
  (let*[(zte(/ (+ j0 wte) 2.)) ; the point in z space that maps to the trailing edge in w space
        (a (- zte c0))         ; the radius of the circle that will generate our airfoil (the complex bit introduces some rotation)
        (b (- zte j0))         ; analogous to a, mapping a circle of radius b about j0 will result in a flat plate.
        (b^2 (* b b))]         ; We never actually use b directly, and its sign is ambiguous when you derive ther relations, anyway.
    (values zte a b^2)))



; a general form of the Joukowski transform
(define (JofZ j0 b^2)
  (λ (z) (+ z (/ b^2 (- z j0)))))

; and its derivative
(define (dJ/dZ j0 b^2)
  (λ (z) (- 1. (/ b^2 z z))))
  
(module+ test
  (printf "joukowsky.rkt: Tests running.~n")
  (require plot)
  ; free parameters
  (define c0 49.+3.i)
  (define j0 50.+1.i)
  (define wte 100.)
  ; dependent parameters
  (define-values (zte a b^2) (make-dependent-parameters c0 j0 wte))
  ; functions to be plotted
  (define c (make-circle c0 a))
  (define t (JofZ j0 b^2))
  ; plotting helper
  (define (z->xy z)(vector(real-part z)(imag-part z)))
  (plot (list (axes)
              (points (list (z->xy c0)(z->xy j0)(z->xy wte)))
              (parametric (λ (theta) (z->xy (c theta))) 0. (* 2. pi))
              (parametric (λ (theta) (z->xy (t (c theta)))) 0. (* 2. pi))))
  (printf "A circle and a reasonable looking airfoil (both stretched vertically) should be displayed.~n"))

;(module+ main
;  (require (submod ".." test)))