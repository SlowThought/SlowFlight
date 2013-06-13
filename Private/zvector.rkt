#lang racket
;;; zvector.rkt provides some common vector operations using complex numbers
;;; to represent the vectors
(provide [contract-out (zdot (-> complex? complex? real?))
                       (zcross (-> complex? complex? real?))])

(require "test-private.rkt")

; dot product
(define (zdot z1 z2)
  (+(*(real-part z1)(real-part z2))
    (*(imag-part z1)(imag-part z2))))

(module+ test
  (printf "zvector running tests.~n")
  (check-eq? (zdot 1 0+i) 0)
  (check-eq? (zdot 0+i 0-i) -1))

; cross product
(define (zcross z1 z2)
  (-(*(real-part z1)(imag-part z2))
    (*(imag-part z1)(real-part z2))))

(module+ test
  (check-eq? (zcross 1+i -1-i) 0)
  (check-eq? (zcross 2 0+3i) 6))