#lang racket/base
(provide integrate)
(define(integrate f0 f1 n)
  (let[(n+1(add1 n))]
    (/(-(expt f1 n+1)
        (expt f0 n+1))
      n+1)))
(module+ test
  (require "test-private.rkt")
  ; No reason to assume that f isn't complex - might come in handy
  ; int(z) = 1/2 z^2. Eval between i and 1.
  (check-complex-= (integrate 0+1.i 1. 1) 1. epsilon))