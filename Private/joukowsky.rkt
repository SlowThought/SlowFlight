#lang racket/base
(provide (all-defined-out))
(require racket/math)
#| The Joukowsky transform maps circles near the origin into airfoil-like shapes, while
   leaving the far field undistorted.

   The parameter a is roughly equivalent to 1/4 of the chord of the airfoil of interest.
   c is a point on the circle C to be transformed. Note that C must contain a. |#
(define(J a c)(+ c(/(sqr a)c)))
(define(dJ/dc a c)(- 1(/(sqr a)c c)))
(module+ test
  (printf "joukowsky.rkt: Tests running.~n")
  (require "test-private.rkt")
  (for-each(λ(c n d)
             (check-=(J 1 c)n epsilon)
             (check-=(dJ/dc 1 c)d epsilon))
           '(1 0+i -1 0-i)
           '(2 0  -2 0)
           '(0 2 0 2)))
(module+ main
  (require (submod ".." test)))