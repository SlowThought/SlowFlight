#lang racket/base
(provide integrate)
#| f is a linear function between x1 & x2 (f = mx + b).
   Find definite integral of f^n in the interval [x1, x2]. |#

(define(integrate f1 f2 x1 x2 n)
  (define dx (- x2 x1))
  (if (= dx 0.) 0.
      (let [(m (/ (- f2 f1) dx))]
        (if (= m 0) 
            (* (expt f1 n) dx)
            (let [(n+1 (add1 n))]
              (/ (- (expt f2 n+1)(expt f1 n+1))
                 m n+1))))))
 
(module+ test
  (require "test-private.rkt")
  (printf "quad.rkt running tests.~n")
  ; delta x = 0
  (check-= (integrate 1. 2. 3. 3. 7) 0. epsilon)
  ; f = 1, F = x
  (check-= (integrate 1. 1. 2. 3. 1) 1. epsilon)
  (check-= (integrate 2. 2. 0. 1. 2) 4. epsilon)
  (check-= (integrate 3. 3. 6. 7. 3) 27. epsilon)
  ; f = x + b, F = 1/2 x^2 + bx
  (check-= (integrate 0. 1. 0. 1. 1) .5 epsilon)
  (check-= (integrate 0. 2. 0. 2. 1) 2. epsilon)
  (check-= (integrate 0. 2. 2. 4. 1) 2. epsilon)
  (check-= (integrate 0. 1. 0. 1. 2) (/ 3.) epsilon)
  (check-= (integrate 0. 1. 7. 8. 2) (/ 3.) epsilon))